// Fixed local validation/numeric regressions. The 77 source-reference cases
// are owned by parity-tests/tests/special_regression_structure_tags.rs.
use std::process::Command;

use cosmolkit_core::{
    StereoError, StructureTagParams, ValenceAssignment, assign_chiral_tags_from_structure,
};
use cosmolkit_model::PropertyValue;
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer3D, CoordinateBlock,
    TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, ChiralTag, Element};

const NUMERIC_CHILD_ENV: &str = "COSMOLKIT_STRUCTURE_TAG_NUMERIC_CHILD";
const NONTETRAHEDRAL_ENV: &str = "RDK_ENABLE_NONTETRAHEDRAL_STEREO";

fn simple_star(center_atomic_number: u8, explicit_hydrogens: u8) -> TopologyBlock {
    let center = AtomSpec::new(
        Element::from_atomic_number(center_atomic_number).expect("modeled center element"),
    )
    .with_explicit_hydrogens(explicit_hydrogens);
    let atoms = [
        center,
        AtomSpec::new(Element::F),
        AtomSpec::new(Element::CL),
        AtomSpec::new(Element::BR),
    ]
    .into_iter()
    .enumerate()
    .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
    .collect();
    let bonds = (1..4)
        .map(|neighbor| {
            Bond::from_spec(
                BondId::new(neighbor - 1),
                BondSpec::new(AtomId::new(0), AtomId::new(neighbor), BondOrder::Single),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn simple_coordinates(rows: Vec<[f64; 3]>) -> CoordinateBlock {
    CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(0, rows, true)],
        ..CoordinateBlock::default()
    }
}

#[test]
fn model_boundaries_return_precise_topology_valence_and_coordinate_errors() {
    let invalid_topology = TopologyBlock {
        atoms: vec![Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C))],
        adjacency: AdjacencyList::from_topology(1, &[]),
        ..TopologyBlock::default()
    };
    assert_eq!(
        assign_chiral_tags_from_structure(
            &invalid_topology,
            &CoordinateBlock::default(),
            &ValenceAssignment {
                explicit_valence: vec![0],
                implicit_hydrogens: vec![0],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::InvalidTopology(
            TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            }
        ))
    );

    let topology = simple_star(6, 1);
    let coordinates = simple_coordinates(vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]);
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 3],
                implicit_hydrogens: vec![0; 4],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::InvalidValence {
            field: "explicit_valence",
            actual: 3,
            atom_count: 4,
        })
    );
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 4],
                implicit_hydrogens: vec![0; 2],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::InvalidValence {
            field: "implicit_hydrogens",
            actual: 2,
            atom_count: 4,
        })
    );
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0, 0, -1, 0],
                implicit_hydrogens: vec![0; 4],
            },
            &StructureTagParams::default(),
        ),
        // Native assignChiralTypesFrom3D reads getTotalNumHs only; this
        // original negative explicit cache is never read. Require the entire
        // output to equal the same original graph with valid explicit rows.
        assign_chiral_tags_from_structure(
            &topology,
            &coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 4],
                implicit_hydrogens: vec![0; 4]
            },
            &StructureTagParams::default(),
        )
    );
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 4],
                implicit_hydrogens: vec![0, -2, 0, 0],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::InvalidValenceValue {
            field: "implicit_hydrogens",
            atom: AtomId::new(1),
            value: -2,
        })
    );

    let short_coordinates = simple_coordinates(vec![[0.0, 0.0, 0.0]; 3]);
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &short_coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 4],
                implicit_hydrogens: vec![0; 4],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::ConformerAtomCountMismatch {
            conformer: 0,
            rows: 3,
            atom_count: 4,
        })
    );

    let duplicate_coordinates = CoordinateBlock {
        conformers_3d: vec![
            Conformer3D::new(9, vec![[0.0, 0.0, 0.0]; 4], true),
            Conformer3D::new(9, vec![[1.0, 0.0, 0.0]; 4], true),
        ],
        ..CoordinateBlock::default()
    };
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &duplicate_coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 4],
                implicit_hydrogens: vec![0; 4],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::DuplicateConformerId { id: 9 })
    );
}

#[test]
fn wedge_and_dash_explicitness_prevent_nonexplicit_property_insertion() {
    for direction in [BondDirection::BeginWedge, BondDirection::BeginDash] {
        let mut topology = simple_star(6, 1);
        topology.bonds[0].set_direction(direction);
        let before = topology.clone();
        let coordinates = simple_coordinates(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        let assignment = assign_chiral_tags_from_structure(
            &topology,
            &coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 4],
                implicit_hydrogens: vec![0; 4],
            },
            &StructureTagParams::default(),
        )
        .unwrap();
        assert_eq!(
            assignment.topology.atoms[0].chiral_tag(),
            ChiralTag::TetrahedralCcw
        );
        assert_eq!(
            assignment.topology.atoms[0].prop("_NonExplicit3DChirality"),
            None
        );
        assert_eq!(topology, before);
    }
}

#[test]
fn numeric_extremes_run_in_a_process_with_source_default_environment() {
    let output = Command::new(std::env::current_exe().expect("current test executable"))
        .arg("--exact")
        .arg("structure_tags_numeric_extremes_child")
        .arg("--nocapture")
        .env(NUMERIC_CHILD_ENV, "1")
        .env_remove(NONTETRAHEDRAL_ENV)
        .output()
        .expect("run isolated numeric boundary test");
    assert!(
        output.status.success(),
        "numeric boundary child failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn structure_tags_numeric_extremes_child() {
    if std::env::var_os(NUMERIC_CHILD_ENV).is_none() {
        return;
    }
    // Native one-opposite-pair T-shape: pair[0]==2 produces UInt permutation 2.
    let platinum = simple_star(78, 0);
    let t_shape = simple_coordinates(vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ]);
    let result = assign_chiral_tags_from_structure(
        &platinum,
        &t_shape,
        &ValenceAssignment {
            explicit_valence: vec![0; 4],
            implicit_hydrogens: vec![0; 4],
        },
        &StructureTagParams::default(),
    )
    .unwrap();
    assert_eq!(
        result.topology.atoms[0].chiral_tag(),
        ChiralTag::SquarePlanar
    );
    assert_eq!(result.topology.atoms[0].chiral_permutation(), Some(2));
    assert_eq!(
        result.topology.atoms[0].prop("_chiralPermutation"),
        Some(&PropertyValue::UInt(2))
    );
    assert_eq!(
        result.topology.atoms[0].prop("_NonExplicit3DChirality"),
        Some(&PropertyValue::Int(1))
    );
    let tetrahedral = simple_star(6, 1);
    let signed_zero = simple_coordinates(vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -0.0],
    ]);
    let valence = ValenceAssignment {
        explicit_valence: vec![0; 4],
        implicit_hydrogens: vec![0; 4],
    };
    let assignment = assign_chiral_tags_from_structure(
        &tetrahedral,
        &signed_zero,
        &valence,
        &StructureTagParams::default(),
    )
    .unwrap();
    assert_eq!(
        assignment.topology.atoms[0].chiral_tag(),
        ChiralTag::Unspecified
    );

    let phosphorus = simple_star(15, 0);
    let minimum_subnormal = simple_coordinates(vec![
        [0.0, 0.0, 0.0],
        [f64::from_bits(1), 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ]);
    assert_eq!(
        assign_chiral_tags_from_structure(
            &phosphorus,
            &minimum_subnormal,
            &valence,
            &StructureTagParams::default(),
        ),
        Err(StereoError::ZeroLengthVector {
            center: AtomId::new(0),
            neighbor: AtomId::new(1),
        })
    );

    let maximum_finite = simple_coordinates(vec![
        [0.0, 0.0, 0.0],
        [f64::MAX, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ]);
    let assignment = assign_chiral_tags_from_structure(
        &phosphorus,
        &maximum_finite,
        &valence,
        &StructureTagParams::default(),
    )
    .unwrap();
    assert_eq!(
        assignment.topology.atoms[0].chiral_tag(),
        ChiralTag::Unspecified
    );
}

#[test]
fn no_implicit_short_circuit_accepts_uninitialized_cache_and_emits_native_int_marker() {
    let mut topology = simple_star(6, 1);
    for atom in &mut topology.atoms {
        atom.set_no_implicit(true);
    }
    let before = topology.clone();
    let coordinates = simple_coordinates(vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]);
    let assignment = assign_chiral_tags_from_structure(
        &topology,
        &coordinates,
        &ValenceAssignment {
            explicit_valence: vec![-1; 4],
            implicit_hydrogens: vec![-1; 4],
        },
        &StructureTagParams::default(),
    )
    .unwrap();
    assert_eq!(
        assignment.topology.atoms[0].chiral_tag(),
        ChiralTag::TetrahedralCcw
    );
    assert_eq!(
        assignment.topology.atoms[0].prop("_NonExplicit3DChirality"),
        Some(&PropertyValue::Int(1))
    );
    assert!(assignment.clear_stereochem_done);
    assert_eq!(topology, before);
}

#[test]
fn existing_tags_skip_implicit_cache_reads_before_degree_checks() {
    let mut topology = simple_star(6, 1);
    for atom in &mut topology.atoms {
        atom.set_chiral_tag(ChiralTag::TetrahedralCw);
    }
    let before = topology.clone();
    let coordinates = simple_coordinates(vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]);
    let assignment = assign_chiral_tags_from_structure(
        &topology,
        &coordinates,
        &ValenceAssignment {
            explicit_valence: vec![-1; 4],
            implicit_hydrogens: vec![-1; 4],
        },
        &StructureTagParams {
            conformer_id: 0,
            replace_existing_tags: false,
        },
    )
    .unwrap();
    assert_eq!(assignment.topology, before);
    assert!(assignment.clear_stereochem_done);
    assert_eq!(topology, before);
    // The first untagged leaf reads its cache before the native degree<3 skip.
    topology.atoms[1].set_chiral_tag(ChiralTag::Unspecified);
    let before_error = topology.clone();
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &coordinates,
            &ValenceAssignment {
                explicit_valence: vec![-1; 4],
                implicit_hydrogens: vec![-1; 4]
            },
            &StructureTagParams {
                conformer_id: 0,
                replace_existing_tags: false
            }
        ),
        Err(StereoError::InvalidValenceValue {
            field: "implicit_hydrogens",
            atom: AtomId::new(1),
            value: -1
        })
    );
    assert_eq!(topology, before_error);
}

#[test]
fn false_is3d_xyz_row_returns_before_cache_and_coordinate_shape_reads() {
    let topology = simple_star(6, 1);
    let coordinates = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(0, vec![[1.0, 2.0, 3.0]], false)],
        ..Default::default()
    };
    let before = coordinates.clone();
    let assignment = assign_chiral_tags_from_structure(
        &topology,
        &coordinates,
        &ValenceAssignment {
            explicit_valence: vec![],
            implicit_hydrogens: vec![],
        },
        &StructureTagParams::default(),
    )
    .unwrap();
    assert_eq!(assignment.topology, topology);
    assert_eq!(assignment.selected_conformer_id, Some(0));
    assert!(!assignment.clear_stereochem_done);
    assert_eq!(coordinates, before);
}
