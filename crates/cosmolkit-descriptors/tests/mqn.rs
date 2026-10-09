//! MQN fixed source-derived full42 test proposals for p1/ROOT review.
//! Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8, MQN.cpp (BSD).
//! No oracles/generators/upstream reads occur in these tests.

use cosmolkit_core::{RingFindType, RingInfo, ValenceAssignment, ValenceError};
use cosmolkit_descriptors::{
    DescriptorError, DescriptorInput, DescriptorSearchCause, MQN_VERSION, mqns,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, CoordinateBlock,
    Element, MoleculeProperties, TopologyBlock,
};
use cosmolkit_search::QueryMatchContextError;

#[derive(Clone, Debug, PartialEq)]
struct Fixture {
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    valence: ValenceAssignment,
    rings: RingInfo,
}

impl Fixture {
    fn input(&self) -> DescriptorInput<'_> {
        DescriptorInput::new(
            &self.topology,
            &self.coordinates,
            &self.properties,
            &self.valence,
            &self.rings,
        )
    }
}

// Deliberate prepared-read fixtures: initialized zero cache rows and
// no_implicit=true, not a substitute chemistry or valence implementation.
// Scalar mutations below explicitly test supplied cached/state reads.
fn graph(symbols: &[&str], edges: &[(usize, usize, BondOrder, bool)]) -> Fixture {
    let atoms = symbols
        .iter()
        .enumerate()
        .map(|(i, symbol)| {
            Atom::from_spec(
                AtomId::new(i),
                AtomSpec::new(Element::from_symbol(symbol).unwrap()).with_no_implicit(true),
            )
        })
        .collect::<Vec<_>>();
    let bonds = edges
        .iter()
        .enumerate()
        .map(|(i, &(begin, end, order, aromatic))| {
            Bond::from_spec(
                BondId::new(i),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), order).with_aromatic(aromatic),
            )
        })
        .collect::<Vec<_>>();
    let adjacency = AdjacencyList::try_from_topology(atoms.len(), &bonds).unwrap();
    let topology = TopologyBlock {
        atoms,
        bonds,
        adjacency,
        ..Default::default()
    };
    topology.validate().unwrap();
    let n = symbols.len();
    Fixture {
        topology,
        coordinates: CoordinateBlock::default(),
        properties: MoleculeProperties::default().with_name("MQN borrowed-state sentinel"),
        valence: ValenceAssignment {
            explicit_valence: vec![0; n],
            implicit_hydrogens: vec![0; n],
        },
        rings: RingInfo::new(RingFindType::OtherOrUnknown, n, edges.len()),
    }
}

fn star(center: &str, leaf: &str, degree: usize) -> Fixture {
    let mut symbols = vec![leaf; degree + 1];
    symbols[0] = center;
    let edges = (1..=degree)
        .map(|i| (0, i, BondOrder::Single, false))
        .collect::<Vec<_>>();
    graph(&symbols, &edges)
}

fn cycle(size: usize, order: BondOrder, flags: usize) -> Fixture {
    let edges = (0..size)
        .map(|i| (i, (i + 1) % size, order, i < flags))
        .collect::<Vec<_>>();
    let mut f = graph(&vec!["C"; size], &edges);
    let row = (0..size).collect::<Vec<_>>();
    f.rings.add_ring(&row, &row).unwrap();
    f
}

fn expected(name: &str) -> [u32; 42] {
    FROZEN
        .iter()
        .find(|(label, _)| *label == name)
        .unwrap_or_else(|| panic!("missing frozen literal {name}"))
        .1
}

fn check(f: &Fixture, name: &str) {
    let before = f.clone();
    let literal = expected(name);
    let first = mqns(&f.input(), false).unwrap_or_else(|error| panic!("{name}/false: {error:?}"));
    let second = mqns(&f.input(), true).unwrap_or_else(|error| panic!("{name}/true: {error:?}"));
    assert_eq!(first.as_slice(), literal.as_slice(), "{name}/false full42");
    assert_eq!(second.as_slice(), literal.as_slice(), "{name}/true full42");
    assert_ne!(
        first.as_ptr(),
        second.as_ptr(),
        "separate owned output vectors"
    );
    assert_eq!(f, &before, "all five borrowed inputs preserved");
    assert_eq!(
        mqns(&f.input(), false).unwrap(),
        first,
        "repeat deterministic"
    );
    assert_eq!(f, &before);
}

fn error(f: &Fixture) -> DescriptorError {
    let before = f.clone();
    let first = mqns(&f.input(), false).unwrap_err();
    assert_eq!(mqns(&f.input(), true).unwrap_err(), first);
    assert_eq!(f, &before, "failure preserves complete inputs");
    first
}

fn real(smiles: &str) -> Fixture {
    let params = cosmolkit_smiles::SmilesParseParams {
        remove_hs: false,
        ..Default::default()
    };
    let record = cosmolkit_smiles::parse_smiles(smiles, &params).unwrap();
    let valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &record.topology,
        cosmolkit_core::ValenceModel::RdkitLike,
        false,
    )
    .unwrap();
    let rings = cosmolkit_core::symmetrized_sssr(&record.topology, &Default::default()).unwrap();
    Fixture {
        topology: record.topology,
        coordinates: CoordinateBlock::default(),
        properties: MoleculeProperties::default().with_name("real prepared MQN fixture"),
        valence,
        rings,
    }
}

#[test]
fn mqn_version_empty_and_ignored_force() {
    assert_eq!(MQN_VERSION, "1.0.0");
    check(&graph(&[], &[]), "empty");
}

#[test]
fn mqn_element_bins_and_heavy_exclusions() {
    let mut f = graph(
        &[
            "C", "F", "Cl", "Br", "I", "S", "P", "N", "O", "Si", "*", "H", "H",
        ],
        &[],
    );
    f.topology.atoms[12] = Atom::from_spec(
        AtomId::new(12),
        AtomSpec::new(Element::H)
            .with_isotope(2)
            .with_no_implicit(true),
    );
    check(&f, "isolated_elements");
    let mut ring = graph(
        &["C", "N", "O"],
        &[
            (0, 1, BondOrder::Single, false),
            (1, 2, BondOrder::Single, false),
            (2, 0, BondOrder::Single, false),
        ],
    );
    ring.rings.add_ring(&[0, 1, 2], &[0, 1, 2]).unwrap();
    check(&ring, "cyclic_C_N_O");
}

#[test]
fn mqn_nitrogen_degree_not_four_and_donors() {
    for d in 0..=5 {
        let mut f = star("N", "C", d);
        f.topology.atoms[0].set_explicit_hydrogens(2);
        check(&f, &format!("nitrogen_degree_{d}"));
    }
}

#[test]
fn mqn_oxygen_charge_and_neighbor_hydrogens() {
    for charge in [-2, -1, 0, 1, 2] {
        let mut f = graph(&["O"], &[]);
        f.topology.atoms[0] = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::O)
                .with_formal_charge(charge)
                .with_explicit_hydrogens(3)
                .with_no_implicit(true),
        );
        check(&f, &format!("oxygen_charge_{charge}"));
    }
    // getTotalNumHs(false) omits explicit H-neighbor rows, even isotope H.
    let mut f = star("O", "H", 2);
    f.topology.atoms[2] = Atom::from_spec(
        AtomId::new(2),
        AtomSpec::new(Element::H)
            .with_isotope(2)
            .with_no_implicit(true),
    );
    check(&f, "oxygen_two_explicit_H_neighbors");
}

#[test]
fn mqn_charge_signs_include_hydrogen_and_wildcard() {
    let mut f = graph(&["*", "H", "Si"], &[]);
    for (i, (element, charge)) in [(Element::DUMMY, -3), (Element::H, 2), (Element::SI, -1)]
        .into_iter()
        .enumerate()
    {
        f.topology.atoms[i] = Atom::from_spec(
            AtomId::new(i),
            AtomSpec::new(element)
                .with_formal_charge(charge)
                .with_no_implicit(true),
        );
    }
    check(&f, "charged_wildcard_H_Si");
}

#[test]
fn mqn_explicit_degrees_and_ring_splits() {
    for d in 1..=5 {
        check(&star("C", "*", d), &format!("carbon_degree_{d}"));
    }
    for d in 2..=4 {
        let symbols = vec!["C"; d + 1];
        let mut edges = vec![
            (0, 1, BondOrder::Single, false),
            (1, 2, BondOrder::Single, false),
            (2, 0, BondOrder::Single, false),
        ];
        for i in 3..=d {
            edges.push((0, i, BondOrder::Single, false));
        }
        let mut f = graph(&symbols, &edges);
        f.rings.add_ring(&[0, 1, 2], &[0, 1, 2]).unwrap();
        check(&f, &format!("cyclic_degree_{d}"));
    }
}

#[test]
fn mqn_all_twenty_two_bond_orders() {
    use BondOrder::*;
    for order in [
        Unspecified,
        Single,
        Double,
        Triple,
        Quadruple,
        Quintuple,
        Hextuple,
        OneAndHalf,
        TwoAndHalf,
        ThreeAndHalf,
        FourAndHalf,
        FiveAndHalf,
        Aromatic,
        Ionic,
        Hydrogen,
        ThreeCenter,
        DativeOne,
        Dative,
        DativeLeft,
        DativeRight,
        Other,
        Zero,
    ] {
        let name = match order {
            Single => "bond_single",
            Double => "bond_double",
            Triple => "bond_triple",
            _ => "bond_default_order",
        };
        check(&graph(&["*", "*"], &[(0, 1, order, false)]), name);
    }
}

#[test]
fn mqn_single_double_triple_ring_and_acyclic() {
    let edges = [
        (0, 1, BondOrder::Single, false),
        (1, 2, BondOrder::Double, false),
        (2, 3, BondOrder::Triple, false),
    ];
    check(&graph(&["C", "C", "C", "C"], &edges), "acyclic_S_D_T");
    let mut f = graph(
        &["C", "C", "C"],
        &[
            (0, 1, BondOrder::Single, false),
            (1, 2, BondOrder::Double, false),
            (2, 0, BondOrder::Triple, false),
        ],
    );
    f.rings.add_ring(&[0, 1, 2], &[0, 1, 2]).unwrap();
    check(&f, "cyclic_S_D_T");
}

#[test]
fn mqn_aromatic_flags_remainder_and_typed_orders() {
    for count in 0..=3 {
        check(
            &cycle(3, BondOrder::Aromatic, count),
            &format!("aromatic_triangle_flags_{count}"),
        );
    }
    for order in [BondOrder::Single, BondOrder::Double, BondOrder::Triple] {
        check(&cycle(3, order, 3), &format!("aromatic_triangle_{order:?}"));
    }
    for order in [
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Triple,
        BondOrder::Aromatic,
        BondOrder::Other,
    ] {
        // Aromatic FLAG is counted even without supplied ring membership;
        // typed S/D/T bins are retained alongside the separate flag split.
        check(
            &graph(&["*", "*"], &[(0, 1, order, true)]),
            &format!("flagged_acyclic_{order:?}"),
        );
    }
}

#[test]
fn mqn_ring_sizes_three_to_nine_and_ten_or_more() {
    for size in 3..=11 {
        check(
            &cycle(size, BondOrder::Single, 0),
            &format!("cycle_size_{size}"),
        );
    }
}

#[test]
fn mqn_multiply_fused_atoms_and_bonds() {
    let mut f = graph(
        &["C", "C", "C", "C"],
        &[
            (0, 1, BondOrder::Single, false),
            (1, 2, BondOrder::Single, false),
            (2, 0, BondOrder::Single, false),
            (1, 3, BondOrder::Single, false),
            (3, 0, BondOrder::Single, false),
        ],
    );
    f.rings.add_ring(&[0, 1, 2], &[0, 1, 2]).unwrap();
    f.rings.add_ring(&[0, 1, 3], &[0, 3, 4]).unwrap();
    check(&f, "two_fused_triangles");
}

#[test]
fn mqn_supplied_empty_duplicate_and_small_ring_rows() {
    let mut f = cycle(4, BondOrder::Single, 0);
    f.rings = RingInfo::new(RingFindType::OtherOrUnknown, 4, 4);
    check(&f, "square_empty_supplied_rings");
    f.rings.add_ring(&[0, 1, 2, 3], &[0, 1, 2, 3]).unwrap();
    f.rings.add_ring(&[0, 1, 2, 3], &[0, 1, 2, 3]).unwrap();
    check(&f, "square_duplicate_supplied_rings");
    let mut empty = graph(&[], &[]);
    empty.rings.add_ring(&[], &[]).unwrap();
    check(&empty, "empty_ring_row_on_empty_graph");
    let mut small = graph(&["C", "C"], &[(0, 1, BondOrder::Single, false)]);
    small.rings.add_ring(&[0], &[0]).unwrap();
    check(&small, "size_one_ring_row");
    small.rings = RingInfo::new(RingFindType::OtherOrUnknown, 2, 1);
    small.rings.add_ring(&[0, 1], &[0, 0]).unwrap();
    check(&small, "size_two_ring_row");
}

#[test]
fn mqn_real_prepared_default_rotatable_chain() {
    check(&real("CCCC"), "real_CCCC");
}

#[test]
fn mqn_pinned_github623_full_vectors() {
    check(&real("CC*"), "github623_CCstar");
    check(&real("[2H][2H]"), "github623_deuterium_pair");
}

#[test]
fn mqn_invalid_prepared_shapes_and_topology() {
    let valid = graph(&["C", "C"], &[(0, 1, BondOrder::Single, false)]);
    for field in ["explicit_valence", "implicit_hydrogens"] {
        let mut f = valid.clone();
        if field == "explicit_valence" {
            f.valence.explicit_valence.pop();
        } else {
            f.valence.implicit_hydrogens.pop();
        }
        assert_eq!(
            error(&f),
            DescriptorError::InvalidValenceRows {
                function: "mqns",
                field,
                actual: 1,
                expected: 2
            }
        );
    }
    let mut f = valid.clone();
    f.topology.atoms[0] = f.topology.atoms[0].clone().with_id(AtomId::new(7));
    let source = f.topology.validate().unwrap_err();
    assert_eq!(
        error(&f),
        DescriptorError::Search {
            function: "mqns",
            source: DescriptorSearchCause::Context(QueryMatchContextError::InvalidTopology(source))
        }
    );
    for (atoms, bonds, field, actual, expected) in [(1, 1, "atoms", 1, 2), (2, 0, "bonds", 0, 1)] {
        f = valid.clone();
        f.rings = RingInfo::new(RingFindType::OtherOrUnknown, atoms, bonds);
        assert_eq!(
            error(&f),
            DescriptorError::Search {
                function: "mqns",
                source: DescriptorSearchCause::Context(
                    QueryMatchContextError::RingMembershipRows {
                        field,
                        actual,
                        expected
                    }
                )
            }
        );
    }
    // Public RingInfo mutators let a caller provide an out-of-range row,
    // then restore dimension lengths. The shared validator retains its cause.
    f = valid.clone();
    f.rings.add_ring(&[7], &[0]).unwrap();
    f.rings.preallocate(2, 1);
    assert!(matches!(
        error(&f),
        DescriptorError::Search {
            function: "mqns",
            source: DescriptorSearchCause::Context(QueryMatchContextError::Rings(
                cosmolkit_core::RingFindingError::RingAtomOutOfRange {
                    atom: 7,
                    atom_count: 2
                }
            ))
        }
    ));
    f = valid;
    f.rings.add_ring(&[0], &[7]).unwrap();
    f.rings.preallocate(2, 1);
    let failure = error(&f);
    assert!(std::error::Error::source(&failure).is_some());
    assert!(matches!(
        failure,
        DescriptorError::Search {
            function: "mqns",
            source: DescriptorSearchCause::Context(QueryMatchContextError::Rings(
                cosmolkit_core::RingFindingError::RingBondOutOfRange {
                    bond: 7,
                    bond_count: 1
                }
            ))
        }
    ));
}

#[test]
fn mqn_prepared_cache_errors_noimplicit_and_preservation() {
    let mut f = graph(&["N"], &[]);
    f.topology.atoms[0].set_explicit_hydrogens(2);
    f.topology.atoms[0].set_no_implicit(false);
    f.valence.implicit_hydrogens[0] = -1;
    let failure = error(&f);
    assert_eq!(
        failure,
        DescriptorError::Valence {
            function: "mqns",
            source: ValenceError::ImplicitValenceCacheNotInitialized {
                atom: AtomId::new(0)
            }
        }
    );
    assert!(
        std::error::Error::source(&failure)
            .unwrap()
            .downcast_ref::<ValenceError>()
            .is_some()
    );
    f.topology.atoms[0].set_no_implicit(true);
    check(&f, "nitrogen_noimplicit_explicit_2");
    f.topology.atoms[0].set_no_implicit(false);
    f.topology.atoms[0].set_explicit_hydrogens(1);
    f.valence.implicit_hydrogens[0] = 2;
    check(&f, "nitrogen_supplied_implicit_2_explicit_1");
    let mut large = graph(&["N", "N", "N"], &[]);
    for atom in &mut large.topology.atoms {
        atom.set_no_implicit(false);
    }
    large.valence.implicit_hydrogens.fill(i32::MAX);
    // Preserve the original three-N wide-cache input as the source's signed
    // byte initialization failure, including both force policies and inputs.
    assert_eq!(
        error(&large),
        DescriptorError::Valence {
            function: "mqns",
            source: ValenceError::ImplicitValenceCacheNotInitialized {
                atom: AtomId::new(0)
            }
        }
    );
    // NoImplicit bypasses even that negative stored cache. All 42 bins remain
    // independently literal, and every original frozen vector stays retained.
    for atom in &mut large.topology.atoms {
        atom.set_no_implicit(true);
    }
    check(&large, "three_N_noimplicit_negative_cache");
    for atom in &mut large.topology.atoms {
        atom.set_no_implicit(false);
        atom.set_explicit_hydrogens(255);
    }
    large.valence.implicit_hydrogens.fill(i32::MAX - 128); // stored int8_t 127
    check(&large, "three_N_maximum_cached_hydrogens");
    for bad_atom in 0..3 {
        large.valence.implicit_hydrogens.fill(127);
        large.valence.implicit_hydrogens[bad_atom] = i32::MAX;
        assert_eq!(
            error(&large),
            DescriptorError::Valence {
                function: "mqns",
                source: ValenceError::ImplicitValenceCacheNotInitialized {
                    atom: AtomId::new(bad_atom)
                }
            }
        );
    }
}

// Full literal vectors frozen in lane receipts before test writing.
const FROZEN: [(&str, [u32; 42]); 65] = [
    (
        "empty",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "isolated_elements",
        [
            1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 10, 0, 0, 0, 0, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "nitrogen_degree_0",
        [
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "nitrogen_degree_1",
        [
            1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "nitrogen_degree_2",
        [
            2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 2, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "nitrogen_degree_3",
        [
            3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 4, 3, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 3, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "nitrogen_degree_4",
        [
            4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 4, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "nitrogen_degree_5",
        [
            5, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 6, 5, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 5, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "oxygen_charge_-2",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 1, 3, 1, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "oxygen_charge_-1",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 3, 1, 3, 1, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "oxygen_charge_0",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 1, 3, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "oxygen_charge_1",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 1, 3, 1, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "oxygen_charge_2",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 1, 3, 1, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "oxygen_two_explicit_H_neighbors",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "charged_wildcard_H_Si",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "carbon_degree_1",
        [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "carbon_degree_2",
        [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "carbon_degree_3",
        [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "carbon_degree_4",
        [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "carbon_degree_5",
        [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "cyclic_degree_2",
        [
            3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "cyclic_degree_3",
        [
            4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
            2, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "cyclic_degree_4",
        [
            5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 2, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "bond_default_order",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "bond_single",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "bond_double",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "bond_triple",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "acyclic_S_D_T",
        [
            4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "cyclic_S_D_T",
        [
            3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "aromatic_triangle_flags_0",
        [
            3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "aromatic_triangle_flags_1",
        [
            3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "aromatic_triangle_flags_2",
        [
            3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "aromatic_triangle_flags_3",
        [
            3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "aromatic_triangle_Single",
        [
            3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "aromatic_triangle_Double",
        [
            3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "aromatic_triangle_Triple",
        [
            3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 2, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "flagged_acyclic_Single",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "flagged_acyclic_Double",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "flagged_acyclic_Triple",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "flagged_acyclic_Aromatic",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "flagged_acyclic_Other",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "cycle_size_3",
        [
            3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "cycle_size_4",
        [
            4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "cycle_size_5",
        [
            5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            5, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "cycle_size_6",
        [
            6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            6, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "cycle_size_7",
        [
            7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "cycle_size_8",
        [
            8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            8, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        ],
    ),
    (
        "cycle_size_9",
        [
            9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            9, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        ],
    ),
    (
        "cycle_size_10",
        [
            10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        ],
    ),
    (
        "cycle_size_11",
        [
            11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 0, 0, 0, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        ],
    ),
    (
        "two_fused_triangles",
        [
            4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            2, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 1,
        ],
    ),
    (
        "square_empty_supplied_rings",
        [
            4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "square_duplicate_supplied_rings",
        [
            4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            4, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 4,
        ],
    ),
    (
        "empty_ring_row_on_empty_graph",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "size_one_ring_row",
        [
            2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "size_two_ring_row",
        [
            2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        ],
    ),
    (
        "real_CCCC",
        [
            4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 3, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "github623_CCstar",
        [
            2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "github623_deuterium_pair",
        [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "nitrogen_noimplicit_explicit_2",
        [
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "nitrogen_supplied_implicit_2_explicit_1",
        [
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 3, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "cyclic_C_N_O",
        [
            1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 3, 0, 0, 0, 3, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "three_N_donor_unsigned_wrap",
        [
            0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 3, 3, 2147483645, 3, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "three_N_noimplicit_negative_cache",
        [
            0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "three_N_maximum_cached_hydrogens",
        [
            0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 3, 3, 1146, 3, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
];
