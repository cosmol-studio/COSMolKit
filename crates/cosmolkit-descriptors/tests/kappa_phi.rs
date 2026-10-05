//! Proposed fixed regressions for p1 review and ROOT adjudication.
//! Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8 (BSD):
//! ConnectivityDescriptors.cpp:297-361; ROMol.cpp:187-196,299-313;
//! Subgraphs.cpp:190-280,451-555; Descriptors/test.cpp:1366-1474.
//! No reference execution, generators or upstream-file reads occur in tests.

use cosmolkit_core::{PathError, PathSearchParams, all_paths_of_length};
use cosmolkit_descriptors::{
    DescriptorError, DescriptorResult, KAPPA_1_VERSION, KAPPA_2_VERSION, KAPPA_3_VERSION,
    PHI_VERSION, kappa_1, kappa_2, kappa_3, phi,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element,
    Hybridization, TopologyBlock, TopologyValidationError,
};

type Entry = fn(&TopologyBlock) -> DescriptorResult<f64>;
const ENTRIES: [(&str, Entry); 4] = [
    ("kappa_1", kappa_1),
    ("kappa_2", kappa_2),
    ("kappa_3", kappa_3),
    ("phi", phi),
];

fn graph(symbols: &[&str], edges: &[(usize, usize)]) -> TopologyBlock {
    let atoms = symbols
        .iter()
        .enumerate()
        .map(|(i, symbol)| {
            Atom::from_spec(
                AtomId::new(i),
                AtomSpec::new(Element::from_symbol(symbol).unwrap()),
            )
        })
        .collect::<Vec<_>>();
    let bonds = edges
        .iter()
        .enumerate()
        .map(|(i, &(begin, end))| {
            Bond::from_spec(
                BondId::new(i),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
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
    topology
}

fn chain(symbols: &[&str]) -> TopologyBlock {
    let edges = (1..symbols.len()).map(|i| (i - 1, i)).collect::<Vec<_>>();
    graph(symbols, &edges)
}

fn values(topology: &TopologyBlock) -> [f64; 4] {
    ENTRIES.map(|(name, entry)| {
        let before = topology.clone();
        let value = entry(topology).unwrap_or_else(|error| panic!("{name}: {error}"));
        assert_eq!(topology, &before, "{name}: borrowed input preservation");
        value
    })
}

fn assert_close(topology: &TopologyBlock, expected: [f64; 4]) {
    for ((name, _), (actual, expected)) in ENTRIES
        .into_iter()
        .zip(values(topology).into_iter().zip(expected))
    {
        assert!(
            (actual - expected).abs() < 1e-12,
            "{name}: {actual:?} vs {expected:?}"
        );
    }
}

fn assert_bits(topology: &TopologyBlock, expected: [f64; 4]) {
    for ((name, _), (actual, expected)) in ENTRIES
        .into_iter()
        .zip(values(topology).into_iter().zip(expected))
    {
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "{name}: {actual:?} vs {expected:?}"
        );
    }
}

#[test]
fn kappa_phi_pinned_versions() {
    assert_eq!(
        [KAPPA_1_VERSION, KAPPA_2_VERSION, KAPPA_3_VERSION],
        ["1.1.0"; 3]
    );
    assert_eq!(PHI_VERSION, "1.0.0");
}

#[test]
fn kappa_phi_fixed_acyclic_graphs_and_heavy_parity() {
    // Fixed graph observations from the pinned helper equations, not values
    // obtained from either the implementation or an oracle during testing.
    assert_bits(&TopologyBlock::default(), [0.0; 4]);
    assert_bits(&chain(&["C"]), [0.0; 4]);
    assert_bits(&chain(&["C", "C"]), [2.0, 0.0, 0.0, 0.0]);
    assert_bits(&chain(&["C", "C", "C"]), [3.0, 2.0, 0.0, 2.0]);
    assert_bits(&chain(&["C", "C", "C", "C"]), [4.0, 3.0, 2.0, 3.0]);
    assert_bits(&chain(&["C", "C", "C", "C", "C"]), [5.0, 4.0, 4.0, 4.0]);
    assert_bits(
        &chain(&["C", "C", "C", "C", "C", "C"]),
        [6.0, 5.0, 4.0, 5.0],
    );
    assert_bits(&graph(&["C", "C"], &[]), [0.0; 4]);
    let star = graph(&["C", "C", "C", "C"], &[(0, 1), (0, 2), (0, 3)]);
    assert_close(&star, [4.0, 4.0 / 3.0, 0.0, 4.0 / 3.0]);
}

#[test]
fn kappa_phi_cycle_closure_bond_dedup_and_non_shortest_defaults() {
    let triangle = graph(&["C", "C", "C"], &[(0, 1), (1, 2), (2, 0)]);
    let square = graph(&["C", "C", "C", "C"], &[(0, 1), (1, 2), (2, 3), (3, 0)]);
    let params = PathSearchParams::default();
    assert_eq!(all_paths_of_length(&triangle, 2, &params).unwrap().len(), 3);
    // Closing the triangle uses all three bonds; rotations/reversals
    // deduplicate to exactly one bond set rather than zero or six paths.
    assert_eq!(all_paths_of_length(&triangle, 3, &params).unwrap().len(), 1);
    assert_eq!(all_paths_of_length(&square, 2, &params).unwrap().len(), 4);
    assert_eq!(all_paths_of_length(&square, 3, &params).unwrap().len(), 4);
    assert_close(&triangle, [4.0 / 3.0, 2.0 / 9.0, 0.0, 8.0 / 81.0]);
    assert_bits(&square, [2.25, 0.75, 0.125, 0.421875]);
}

#[test]
fn kappa_phi_explicit_and_isotopic_h_bonds_have_distinct_count_roles() {
    let mut topology = chain(&["H", "C", "C", "C", "C"]);
    assert_bits(&topology, [2.25, 3.0, 2.0, 1.6875]);
    topology.atoms[0] = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::H).with_isotope(2));
    assert_bits(&topology, [2.25, 3.0, 2.0, 1.6875]);
    assert_close(
        &chain(&["H", "C", "C", "C"]),
        [4.0 / 3.0, 2.0, 0.0, 8.0 / 9.0],
    );
    assert_bits(&chain(&["H", "H", "H", "H"]), [0.0; 4]);
    // Atom hydrogen properties create no extra explicit P1 bonds.
    let mut decorated = chain(&["C", "C", "C", "C"]);
    decorated.atoms[0] = Atom::from_spec(
        AtomId::new(0),
        AtomSpec::new(Element::C).with_explicit_hydrogens(3),
    );
    assert_bits(&decorated, [4.0, 3.0, 2.0, 3.0]);
}

#[test]
fn kappa_phi_wildcards_enter_paths_but_not_heavy_count() {
    assert_close(&chain(&["C", "*", "*", "C"]), [2.0 / 9.0, 0.0, 0.0, 0.0]);
    assert_bits(&chain(&["*", "C", "C"]), [0.5, 0.0, 0.0, 0.0]);
    assert_bits(&chain(&["*", "*", "*"]), [0.0, -4.0, 0.0, 0.0]);
    // Only Phi has the zero-heavy early return. Kappa2/3 retain negative
    // source results from nonempty wildcard paths; no blanket guard/clamp.
    assert_bits(&chain(&["*", "*", "*", "*"]), [0.0, -1.0, -18.0, 0.0]);
    assert_bits(&graph(&["*", "H", "*"], &[]), [0.0; 4]);
}

#[test]
fn kappa_phi_alpha_denominator_cancellation_negative_and_signed_zero() {
    // Legitimate source rb0=0 for Cm gives alpha=-1; no radius fallback.
    assert_bits(&chain(&["Cm"]), [0.0, -4.0, -9.0, -0.0]);
    assert_bits(&chain(&["C", "Cm"]), [0.0, 0.0, -4.0, 0.0]);
    assert_bits(&chain(&["C", "C", "Cm"]), [2.0, 0.0, 1.0, 0.0]);
    assert_bits(&chain(&["C", "C", "C", "Cm"]), [3.0, 2.0, 0.0, 1.5]);
    // rb0(Na)=1.54, carbon rb0=0.77, alpha=+1. Missing paths do not
    // automatically zero the Kappa3 result when the denominator is nonzero.
    assert_bits(&chain(&["Na"]), [2.0, 0.0, 1.0, 0.0]);
    assert_bits(&chain(&["C", "C", "C", "Na"]), [5.0, 4.0, 3.0, 5.0]);
}

#[test]
fn kappa_phi_stored_hybridization_and_unread_metadata() {
    let mut topology = chain(&["C", "C", "C", "C"]);
    topology.atoms[1] = Atom::from_spec(
        AtomId::new(1),
        AtomSpec::new(Element::C)
            .with_isotope(13)
            .with_formal_charge(1)
            .with_explicit_hydrogens(2)
            .with_aromatic(true),
    );
    // Stored UNSPECIFIED stays the source default even with these metadata.
    assert_bits(&topology, [4.0, 3.0, 2.0, 3.0]);
    topology.atoms[1].set_hybridization(Hybridization::Sp2);
    assert_close(&topology, [3.87, 2.87, 1.87, 2.776725]);
    topology.atoms[1].set_hybridization(Hybridization::Sp);
    assert_close(&topology, [3.78, 2.78, 1.78, 2.6271]);
}

#[test]
fn kappa_phi_malformed_topology_retains_typed_causes_and_input() {
    let valid = chain(&["C", "C"]);
    let mut malformed = Vec::new();
    let mut atom_id = valid.clone();
    atom_id.atoms[0] = atom_id.atoms[0].clone().with_id(AtomId::new(7));
    malformed.push(atom_id);
    let mut bond_id = valid.clone();
    bond_id.bonds[0] = Bond::from_spec(
        BondId::new(7),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
    );
    malformed.push(bond_id);
    let mut endpoint = valid.clone();
    endpoint.bonds[0] = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(9), BondOrder::Single),
    );
    malformed.push(endpoint);
    let mut adjacency = valid.clone();
    adjacency.adjacency = AdjacencyList::default();
    malformed.push(adjacency);
    let mut self_loop = valid.clone();
    self_loop.bonds[0] = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(0), BondOrder::Single),
    );
    malformed.push(self_loop);
    for topology in malformed {
        let expected_source = topology.validate().unwrap_err();
        let before = topology.clone();
        for (name, entry) in ENTRIES {
            let error = entry(&topology).unwrap_err();
            let expected = if name == "kappa_2" || name == "kappa_3" {
                DescriptorError::Path {
                    function: name,
                    source: PathError::InvalidTopology(expected_source.clone()),
                }
            } else {
                DescriptorError::InvalidTopology {
                    function: name,
                    source: expected_source.clone(),
                }
            };
            assert_eq!(error, expected);
            let cause = std::error::Error::source(&error).expect("structural source retained");
            if name == "kappa_2" || name == "kappa_3" {
                assert!(cause.downcast_ref::<PathError>().is_some());
            } else {
                assert!(cause.downcast_ref::<TopologyValidationError>().is_some());
            }
            assert_eq!(topology, before);
        }
    }
}

#[test]
fn kappa_phi_repeated_calls_preserve_all_borrowed_rows() {
    let topology = graph(&["C", "C", "C", "C"], &[(0, 1), (1, 2), (2, 3), (3, 0)]);
    let before = topology.clone();
    for _ in 0..3 {
        assert_bits(&topology, [2.25, 0.75, 0.125, 0.421875]);
        assert_eq!(topology, before);
    }
}

fn source_fixture(smiles: &str, nonstrict: bool) -> TopologyBlock {
    let record = cosmolkit_smiles::parse_smiles(smiles, &Default::default()).unwrap();
    let params = if nonstrict {
        use cosmolkit_core::SanitizeOperations as Ops;
        // Upstream test.cpp uses ALL ^ PROPERTIES after a non-strict property
        // cache update. The detached owner performs that update itself when
        // PROPERTIES is absent. Its raw-bit constructor rejects reserved bits,
        // so spell out every installed named stage except PROPERTIES.
        cosmolkit_core::SanitizeParams {
            operations: Ops::CLEANUP
                | Ops::SYMM_RINGS
                | Ops::KEKULIZE
                | Ops::FIND_RADICALS
                | Ops::SET_AROMATICITY
                | Ops::SET_CONJUGATION
                | Ops::SET_HYBRIDIZATION
                | Ops::CLEANUP_CHIRALITY
                | Ops::ADJUST_HS
                | Ops::CLEANUP_ORGANOMETALLICS
                | Ops::CLEANUP_ATROPISOMERS,
        }
    } else {
        Default::default()
    };
    cosmolkit_core::sanitize_topology(&record.topology, &params)
        .unwrap()
        .topology
}

fn golden(entry: Entry, smiles: &str, expected: f64, nonstrict: bool) {
    let topology = source_fixture(smiles, nonstrict);
    let before = topology.clone();
    let actual = entry(&topology).unwrap();
    // The upstream fixed table rounds to three/four decimals and itself
    // uses 0.002. This is a proposed test tolerance, not an acceptance rule.
    assert!(
        (actual - expected).abs() < 0.002,
        "{smiles}: {actual} vs {expected}"
    );
    assert_eq!(topology, before);
}

#[test]
fn kappa_1_seven_pinned_upstream_molecules() {
    let cases = [
        ("C12CC2C3CC13", 2.344),
        ("C1CCC12CC2", 3.061),
        ("C1CCCCC1", 4.167),
        ("CCCCCC", 6.000),
        ("CCC(C)C1CCC(C)CC1", 9.091),
        ("CC(C)CC1CCC(C)CC1", 9.091),
        ("CC(C)C1CCC(C)CCC1", 9.091),
    ];
    assert_eq!(cases.len(), 7);
    for (smiles, expected) in cases {
        golden(kappa_1, smiles, expected, false);
    }
}

#[test]
fn kappa_2_twenty_three_pinned_upstream_molecules() {
    let cases = [
        ("[C+2](C)(C)(C)(C)(C)C", 0.667),
        ("[C+](C)(C)(C)(C)(CC)", 1.240),
        ("C(C)(C)(C)(CCC)", 2.3444),
        ("CC(C)CCCC", 4.167),
        ("CCCCCCC", 6.000),
        ("CCCCCC", 5.000),
        ("CCCCCCC", 6.000),
        ("C1CCCC1", 1.440),
        ("C1CCCC1C", 1.633),
        ("C1CCCCC1", 2.222),
        ("C1CCCCCC1", 3.061),
        ("CCCCC", 4.000),
        ("CC=CCCC", 4.740),
        ("C1=CN=CN1", 0.884),
        ("c1ccccc1", 1.606),
        ("c1cnccc1", 1.552),
        ("n1ccncc1", 1.500),
        ("CCCCF", 3.930),
        ("CCCCCl", 4.290),
        ("CCCCBr", 4.480),
        ("CCC(C)C1CCC(C)CC1", 4.133),
        ("CC(C)CC1CCC(C)CC1", 4.133),
        ("CC(C)C1CCC(C)CCC1", 4.133),
    ];
    assert_eq!(cases.len(), 23);
    for (smiles, expected) in cases {
        golden(kappa_2, smiles, expected, true);
    }
}

#[test]
fn kappa_3_eight_pinned_upstream_molecules() {
    let cases = [
        ("C[C+](C)(C)(C)C(C)(C)C", 2.000),
        ("CCC(C)C(C)(C)(CC)", 2.380),
        ("CCC(C)CC(C)CC", 4.500),
        ("CC(C)CCC(C)CC", 5.878),
        ("CC(C)CCCC(C)C", 8.000),
        ("CCC(C)C1CCC(C)CC1", 2.500),
        ("CC(C)CC1CCC(C)CC1", 3.265),
        ("CC(C)C1CCC(C)CCC1", 2.844),
    ];
    assert_eq!(cases.len(), 8);
    for (smiles, expected) in cases {
        golden(kappa_3, smiles, expected, true);
    }
}
