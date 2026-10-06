//! Public-boundary ring-transport regressions for the sanitize owner.
//!
//! L01: the `final_rings` transport semantics at the two extreme operation
//! masks — `NONE` reports the SOURCE-uninitialized final ring state (the
//! sanitizeMol entry `clearComputedProps()` resets rings for every mask),
//! and `SYMM_RINGS` on an acyclic/empty graph reports the initialized
//! exact final stage state, including initialized-empty rows.

use cosmolkit_core::RingFindType;
use cosmolkit_core::{SanitizeOperations, SanitizeParams, sanitize_topology};
use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
use cosmolkit_types::{BondOrder, Element};

fn topology_from(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn cc_graph() -> TopologyBlock {
    topology_from(
        vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ],
        vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
        )],
    )
}

fn benzene() -> TopologyBlock {
    topology_from(
        (0..6)
            .map(|id| {
                Atom::from_spec(
                    AtomId::new(id),
                    AtomSpec::new(Element::C).with_aromatic(true),
                )
            })
            .collect(),
        (0..6)
            .map(|id| {
                Bond::from_spec(
                    BondId::new(id),
                    BondSpec::new(
                        AtomId::new(id),
                        AtomId::new((id + 1) % 6),
                        BondOrder::Aromatic,
                    )
                    .with_aromatic(true),
                )
            })
            .collect(),
    )
}

fn cube() -> TopologyBlock {
    let atoms = (0..8)
        .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C)))
        .collect::<Vec<_>>();
    let edges = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];
    let bonds = edges
        .iter()
        .enumerate()
        .map(|(id, (begin, end))| {
            Bond::from_spec(
                BondId::new(id),
                BondSpec::new(AtomId::new(*begin), AtomId::new(*end), BondOrder::Single),
            )
        })
        .collect::<Vec<_>>();
    topology_from(atoms, bonds)
}

fn literal_rows(rings: &cosmolkit_core::RingInfo) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    (
        rings
            .atom_rings()
            .iter()
            .map(|row| row.iter().map(|atom| atom.index()).collect())
            .collect(),
        rings
            .bond_rings()
            .iter()
            .map(|row| row.iter().map(|bond| bond.index()).collect())
            .collect(),
    )
}

// L08: the complete 64-call public sanitizer matrix — four graphs x ALL16
// subsets of {SYMM_RINGS, KEKULIZE, SET_AROMATICITY, CLEANUP_ATROPISOMERS}
// (no other bits). Independent literal dispositions; cube Symm rows are
// compared as NORMALIZED face SETS (set/graph evidence, not an ordered
// finder oracle) with the three-per-atom / two-per-bond membership counts.
// One Some disposition = exactly one owning acquisition event of that type;
// None = zero, so the final find_type/rows are the public stage-acquisition
// census at this boundary.
#[test]
fn sanitize_ring_l08_public_sixty_four_call_matrix() {
    use cosmolkit_core::RingFindType;

    fn normalize(row: &[usize]) -> Vec<usize> {
        // Faces are undirected cycles: identity is the SORTED atom
        // multiset; bond correspondence is compared separately as an
        // undirected edge set below.
        let mut sorted = row.to_vec();
        sorted.sort_unstable();
        sorted
    }

    let empty = topology_from(Vec::new(), Vec::new());
    let cc = cc_graph();
    let benzene = benzene();
    let cube = cube();
    const S: u32 = 1;
    const K: u32 = 2;
    const A: u32 = 4;
    const T: u32 = 8;

    let mut calls = 0usize;
    for mask in 0..16u32 {
        let mut operations = SanitizeOperations::NONE;
        if mask & S != 0 {
            operations = operations | SanitizeOperations::SYMM_RINGS;
        }
        if mask & K != 0 {
            operations = operations | SanitizeOperations::KEKULIZE;
        }
        if mask & A != 0 {
            operations = operations | SanitizeOperations::SET_AROMATICITY;
        }
        if mask & T != 0 {
            operations = operations | SanitizeOperations::CLEANUP_ATROPISOMERS;
        }
        let s = mask & S != 0;
        let k = mask & K != 0;
        let a = mask & A != 0;
        let t = mask & T != 0;
        let _ = t; // T never changes the ring disposition; it must not.
        for (name, graph, atom_count, bond_count, expect) in [
            (
                "empty",
                &empty,
                0usize,
                0usize,
                if s || a { "symm-empty" } else { "none" },
            ),
            ("cc", &cc, 2, 1, if s || a { "symm-empty" } else { "none" }),
            (
                "benzene",
                &benzene,
                6,
                6,
                if s {
                    "symm1"
                } else if k {
                    "sssr1"
                } else if a {
                    "symm1"
                } else {
                    "none"
                },
            ),
            ("cube", &cube, 8, 12, if s || a { "symm6" } else { "none" }),
        ] {
            let graph_snapshot = graph.clone();
            let assignment = sanitize_topology(
                graph,
                &SanitizeParams {
                    operations,
                    ..SanitizeParams::default()
                },
            )
            .unwrap_or_else(|error| panic!("{name}/m{mask}: {error:?}"));
            calls += 1;
            assert_eq!(graph, &graph_snapshot, "{name}/m{mask}: input mutated");
            match expect {
                "none" => {
                    assert!(
                        assignment.final_rings.is_none(),
                        "{name}/m{mask}: expected None, got {:?}",
                        assignment.final_rings.as_ref().map(|r| r.find_type())
                    );
                }
                "symm-empty" => {
                    let rings = assignment.final_rings.as_ref().unwrap();
                    assert_eq!(rings.find_type(), RingFindType::SymmSssr, "{name}/m{mask}");
                    assert!(rings.is_initialized(), "{name}/m{mask}");
                    assert!(rings.is_symm_sssr(), "{name}/m{mask}");
                    assert_eq!(rings.atom_row_count(), atom_count, "{name}/m{mask}");
                    assert_eq!(rings.bond_row_count(), bond_count, "{name}/m{mask}");
                    assert!(rings.atom_rings().is_empty(), "{name}/m{mask}: rows");
                }
                "symm1" => {
                    let rings = assignment.final_rings.as_ref().unwrap();
                    assert_eq!(rings.find_type(), RingFindType::SymmSssr, "{name}/m{mask}");
                    assert_eq!(rings.atom_rings().len(), 1, "{name}/m{mask}");
                    let (atoms, bonds) = literal_rows(rings);
                    assert_eq!(atoms, vec![vec![0, 5, 4, 3, 2, 1]], "{name}/m{mask}");
                    assert_eq!(bonds, vec![vec![5, 4, 3, 2, 1, 0]], "{name}/m{mask}");
                }
                "sssr1" => {
                    let rings = assignment.final_rings.as_ref().unwrap();
                    assert_eq!(rings.find_type(), RingFindType::Sssr, "{name}/m{mask}");
                    assert_eq!(rings.atom_rings().len(), 1, "{name}/m{mask}");
                    let (atoms, bonds) = literal_rows(rings);
                    assert_eq!(atoms, vec![vec![0, 5, 4, 3, 2, 1]], "{name}/m{mask}");
                    assert_eq!(bonds, vec![vec![5, 4, 3, 2, 1, 0]], "{name}/m{mask}");
                }
                _ => {
                    // cube Symm: six normalized faces as a SET.
                    let rings = assignment.final_rings.as_ref().unwrap();
                    assert_eq!(rings.find_type(), RingFindType::SymmSssr, "{name}/m{mask}");
                    assert_eq!(rings.atom_rings().len(), 6, "{name}/m{mask}");
                    let (atoms, bonds) = literal_rows(rings);
                    let expected_faces: Vec<(Vec<usize>, Vec<usize>)> = vec![
                        (vec![0, 1, 2, 3], vec![0, 1, 2, 3]),
                        (vec![4, 5, 6, 7], vec![4, 5, 6, 7]),
                        (vec![0, 1, 4, 5], vec![0, 4, 8, 9]),
                        (vec![1, 2, 5, 6], vec![1, 5, 9, 10]),
                        (vec![2, 3, 6, 7], vec![2, 6, 10, 11]),
                        (vec![0, 3, 4, 7], vec![3, 7, 8, 11]),
                    ];
                    // The frozen evidence is per-face SETS: the sorted atom
                    // multiset and the sorted bond multiset of each face.
                    // (Positional pairing would over-constrain the finder's
                    // traversal order, which the packet does not freeze.)
                    let mut got: Vec<(Vec<usize>, Vec<usize>)> = atoms
                        .iter()
                        .zip(&bonds)
                        .map(|(atom_row, bond_row)| (normalize(atom_row), normalize(bond_row)))
                        .collect();
                    let mut expected = expected_faces.clone();
                    got.sort();
                    expected.sort();
                    assert_eq!(got, expected, "{name}/m{mask}: normalized face set");
                    for atom in 0..8usize {
                        assert_eq!(
                            atoms.iter().filter(|row| row.contains(&atom)).count(),
                            3,
                            "{name}/m{mask}: atom {atom} membership"
                        );
                    }
                    for bond in 0..12usize {
                        assert_eq!(
                            bonds.iter().filter(|row| row.contains(&bond)).count(),
                            2,
                            "{name}/m{mask}: bond {bond} membership"
                        );
                    }
                }
            }
        }
    }
    assert_eq!(calls, 64, "exact census: 16 masks x 4 graphs");
}

// L02 four-call benzene matrix: KEKULIZE x {SYMM_RINGS absent/present} x
// {SET_AROMATICITY absent/present}. Absent SYMM leaves the carrier absent so
// the Kekulize owner follows the source uninitialized-state SSSR branch and
// moves its fresh SSSR rows back (final Sssr1). Present SYMM pre-initializes
// the carrier with SymmSssr rows; the Kekulize stage borrows them (final
// SymmSssr1) and AROMATICITY — which today clones the present carrier —
// neither upgrades SSSR to Symm nor re-finds Symm. The find_type of the
// transported state is the public disposition proof: a second permanent
// discovery by K would surface as Sssr even under SYMM_RINGS.
#[test]
fn sanitize_ring_l02_kekulize_four_call_benzene_matrix() {
    let graph = benzene();
    let before = graph.clone();
    struct Cell {
        operations: SanitizeOperations,
        expect_symm: bool,
        aromaticity: bool,
    }
    let cells = [
        Cell {
            operations: SanitizeOperations::KEKULIZE,
            expect_symm: false,
            aromaticity: false,
        },
        Cell {
            operations: SanitizeOperations::KEKULIZE | SanitizeOperations::SET_AROMATICITY,
            expect_symm: false,
            aromaticity: true,
        },
        Cell {
            operations: SanitizeOperations::SYMM_RINGS | SanitizeOperations::KEKULIZE,
            expect_symm: true,
            aromaticity: false,
        },
        Cell {
            operations: SanitizeOperations::SYMM_RINGS
                | SanitizeOperations::KEKULIZE
                | SanitizeOperations::SET_AROMATICITY,
            expect_symm: true,
            aromaticity: true,
        },
    ];
    let mut calls = 0usize;
    for cell in &cells {
        calls += 1;
        let assignment = sanitize_topology(
            &graph,
            &SanitizeParams {
                operations: cell.operations,
            },
        )
        .unwrap();
        assert_eq!(&graph, &before, "input mutated");
        let rings = assignment
            .final_rings
            .as_ref()
            .expect("kekulized benzene transports rows");
        assert!(rings.is_initialized());
        if cell.expect_symm {
            assert_eq!(rings.find_type(), RingFindType::SymmSssr);
        } else {
            assert_eq!(rings.find_type(), RingFindType::Sssr);
        }
        assert_eq!(rings.atom_rings().len(), 1);
        assert_eq!(rings.bond_rings().len(), 1);
        // One six-member row with paired membership dimensions regardless of
        // the finder; the ordered literal belongs to the owning-stage test.
        assert_eq!(rings.atom_rings()[0].len(), 6);
        assert_eq!(rings.bond_rings()[0].len(), 6);
        assert_eq!(rings.atom_row_count(), 6);
        assert_eq!(rings.bond_row_count(), 6);
        for index in 0..6 {
            assert_eq!(rings.atom_members(AtomId::new(index)), &[0]);
            assert_eq!(rings.bond_members(BondId::new(index)), &[0]);
        }
        // Kekulize produced alternating single/double orders and cleared the
        // aromatic flags. SET_AROMATICITY (source semantics: benzene is
        // aromatic after perception) re-marks the ring bonds with the
        // Aromatic type and flag in the A cells.
        if cell.aromaticity {
            for bond in &assignment.topology.bonds {
                assert_eq!(bond.order(), BondOrder::Aromatic);
                assert!(bond.is_aromatic());
            }
        } else {
            let mut doubles = 0;
            for (index, bond) in assignment.topology.bonds.iter().enumerate() {
                assert!(!bond.is_aromatic());
                // K-RING frozen derivation: sanitize's Kekulize runs
                // canonical=false (iota ranks), start atom 0 steps to
                // atom 1 doubling b0: even indexes are Double.
                let expected = if index % 2 == 0 {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                };
                assert_eq!(bond.order(), expected);
                if bond.order() == BondOrder::Double {
                    doubles += 1;
                }
            }
            assert_eq!(doubles, 3);
        }
    }
    assert_eq!(calls, 4, "exact census");
    // Row order/member vectors are identical across the two Symm cells and
    // the two Sssr cells: the ordered literal is a stable property of the
    // owning finders, and the transported state preserves it.
    let symm_a = sanitize_topology(
        &graph,
        &SanitizeParams {
            operations: SanitizeOperations::SYMM_RINGS | SanitizeOperations::KEKULIZE,
        },
    )
    .unwrap();
    let symm_b = sanitize_topology(
        &graph,
        &SanitizeParams {
            operations: SanitizeOperations::SYMM_RINGS
                | SanitizeOperations::KEKULIZE
                | SanitizeOperations::SET_AROMATICITY,
        },
    )
    .unwrap();
    assert_eq!(
        literal_rows(symm_a.final_rings.as_ref().unwrap()),
        literal_rows(symm_b.final_rings.as_ref().unwrap())
    );
}

// L03 eight-call product: {empty, CC, benzene, cube} x {SET_AROMATICITY,
// SYMM_RINGS|SET_AROMATICITY}. The aromaticity guard borrows initialized
// carrier rows (source: srings = atomRings()) and symmetrizes only when the
// carrier is absent — exactly one fresh symmetrized acquisition per A call
// in both mask arms (the SYMM arm pre-seeds the SAME Symm state; A never
// re-finds). Public proof: transported type SymmSssr with identical ordered
// rows across the two arms; empty/CC transport initialized-empty states.
#[test]
fn sanitize_ring_l03_aromaticity_initialized_only_borrowing() {
    let mut calls = 0usize;
    for name in ["empty", "cc", "benzene", "cube"] {
        let graph = match name {
            "empty" => topology_from(vec![], vec![]),
            "cc" => cc_graph(),
            "benzene" => benzene(),
            _ => cube(),
        };
        let before = graph.clone();
        let expected_rows = match name {
            "empty" | "cc" => 0,
            "benzene" => 1,
            _ => 6,
        };
        let mut results = Vec::new();
        for operations in [
            SanitizeOperations::SET_AROMATICITY,
            SanitizeOperations::SYMM_RINGS | SanitizeOperations::SET_AROMATICITY,
        ] {
            calls += 1;
            let assignment = sanitize_topology(&graph, &SanitizeParams { operations }).unwrap();
            assert_eq!(&graph, &before, "{name}: input mutated");
            let rings = assignment
                .final_rings
                .as_ref()
                .expect("{name}: A transports the initialized final state");
            assert!(rings.is_initialized(), "{name}");
            assert_eq!(rings.find_type(), RingFindType::SymmSssr, "{name}");
            assert_eq!(rings.atom_rings().len(), expected_rows, "{name}: rows");
            assert_eq!(rings.bond_rings().len(), expected_rows, "{name}: rows");
            results.push(literal_rows(rings));
        }
        // Borrowed arm identical to the fresh arm: ordered rows/members are
        // the same state, proving no second discovery replaced them.
        assert_eq!(results[0], results[1], "{name}: arms diverged");
        if expected_rows > 0 {
            let rings_shape = &results[0];
            for (atoms, bonds) in rings_shape.0.iter().zip(&rings_shape.1) {
                assert_eq!(atoms.len(), bonds.len(), "{name}: paired rows");
            }
        }
    }
    assert_eq!(calls, 8, "exact census");
}

#[test]
fn sanitize_ring_l01_none_reports_source_uninitialized_final_state() {
    for graph in [topology_from(vec![], vec![]), cc_graph()] {
        let before = graph.clone();
        let assignment = sanitize_topology(
            &graph,
            &SanitizeParams {
                operations: SanitizeOperations::NONE,
            },
        )
        .unwrap();
        // The sanitizeMol entry clearComputedProps() resets ring info for
        // EVERY mask; NONE executes no ring stage, so the final state is
        // the source-uninitialized one — reported as None, distinct from
        // the Kekulize owner's preserve-input None semantics.
        assert!(assignment.final_rings.is_none());
        // Input preservation: the borrowed graph is unchanged.
        assert_eq!(&graph, &before);
    }
}

#[test]
fn sanitize_ring_l01_symm_rings_on_empty_initializes_symm_empty_state() {
    let graph = topology_from(vec![], vec![]);
    let before = graph.clone();
    let assignment = sanitize_topology(
        &graph,
        &SanitizeParams {
            operations: SanitizeOperations::SYMM_RINGS,
        },
    )
    .unwrap();
    let rings = assignment
        .final_rings
        .as_ref()
        .expect("SYMM_RINGS on an empty graph transports an initialized state");
    // Initialized SymmSssr with zero rows and full membership dimensions.
    assert!(rings.is_initialized());
    assert!(rings.is_symm_sssr());
    assert_eq!(rings.atom_rings().len(), 0);
    assert_eq!(rings.bond_rings().len(), 0);
    assert_eq!(rings.atom_row_count(), 0);
    assert_eq!(rings.bond_row_count(), 0);
    assert_eq!(&graph, &before);
}

#[test]
fn sanitize_ring_l01_symm_rings_on_cc_clears_computed_properties_and_preserves_input() {
    // The entry clearComputedProps() drops atom/bond computed properties
    // even when later stages never recompute them; SYMM_RINGS on an
    // acyclic graph transports the initialized-empty Symm state.
    let graph = cc_graph();
    let before = graph.clone();
    let assignment = sanitize_topology(
        &graph,
        &SanitizeParams {
            operations: SanitizeOperations::SYMM_RINGS,
        },
    )
    .unwrap();
    let rings = assignment.final_rings.as_ref().unwrap();
    assert!(rings.is_initialized());
    assert!(rings.is_symm_sssr());
    assert_eq!(rings.atom_rings().len(), 0);
    assert_eq!(&graph, &before);
    // The output topology keeps the input bond order (no aromaticity or
    // kekulize stage ran under this mask).
    assert_eq!(assignment.topology.bonds[0].order(), BondOrder::Double);
}
