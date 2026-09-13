#![cfg(feature = "op-contracts-strict")]

use cosmolkit_core::{
    __migration_hydrogens::{
        AddHydrogensTopologyResult, AddedHydrogen, AddedHydrogenKind, HydrogenError,
        add_hydrogen_coordinates, add_hydrogens_topology,
    },
    AddHsParams, add_hydrogens_with_params, rdkit_rb0,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomPdbResidueInfo, AtomSpec, Bond, BondId, BondOrder, BondSpec,
    ChiralTag, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension,
    CoordinateValidationError, Element, Hybridization, MappingValidationError, MoleculeProperties,
    TopologyBlock, TopologyMapping,
};

const EPSILON: f64 = 1.0e-10;

fn topology(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
    let atoms = atom_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
        .collect::<Vec<_>>();
    let bonds = bond_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
        .collect::<Vec<_>>();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
}

fn explicit_addition(source: TopologyBlock) -> AddHydrogensTopologyResult {
    add_hydrogens_topology(
        source,
        &AddHsParams {
            explicit_only: true,
            ..Default::default()
        },
    )
    .unwrap()
}

fn isolated_parent(explicit_hydrogens: u8, hybridization: Hybridization) -> TopologyBlock {
    topology(
        vec![
            AtomSpec::new(Element::C)
                .with_explicit_hydrogens(explicit_hydrogens)
                .with_hybridization(hybridization),
        ],
        Vec::new(),
    )
}

fn star_parent(other_neighbors: usize, parent_spec: AtomSpec) -> TopologyBlock {
    let mut atoms = vec![parent_spec.with_explicit_hydrogens(1)];
    atoms.extend((0..other_neighbors).map(|_| AtomSpec::new(Element::C)));
    let bonds = (0..other_neighbors)
        .map(|neighbor| BondSpec::new(AtomId::new(0), AtomId::new(neighbor + 1), BondOrder::Single))
        .collect();
    topology(atoms, bonds)
}

fn one_hydrogen_result() -> AddHydrogensTopologyResult {
    explicit_addition(isolated_parent(1, Hybridization::Sp3))
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < EPSILON,
        "{actual} != {expected}"
    );
}

fn distance_3d(left: [f64; 3], right: [f64; 3]) -> f64 {
    ((left[0] - right[0]).powi(2) + (left[1] - right[1]).powi(2) + (left[2] - right[2]).powi(2))
        .sqrt()
}

#[test]
fn empty_and_no_addition_inputs_are_identity_values() {
    let empty = AddHydrogensTopologyResult {
        topology: TopologyBlock::default(),
        mapping: TopologyMapping::identity(0, 0),
        additions: Vec::new(),
    };
    let output =
        add_hydrogen_coordinates(empty.clone(), CoordinateBlock::default(), true, true).unwrap();
    assert_eq!(output.topology, empty.topology);
    assert_eq!(output.mapping, empty.mapping);
    assert!(output.additions.is_empty());
    assert_eq!(output.coordinates, CoordinateBlock::default());

    let source = topology(vec![AtomSpec::new(Element::NE)], Vec::new());
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(5, vec![[4.0, -2.0]])],
        ..Default::default()
    };
    let result = AddHydrogensTopologyResult {
        topology: source.clone(),
        mapping: TopologyMapping::identity(1, 0),
        additions: Vec::new(),
    };
    let output = add_hydrogen_coordinates(result, coordinates.clone(), true, true).unwrap();
    assert_eq!(output.topology, source);
    assert_eq!(output.coordinates, coordinates);
}

#[test]
fn disabled_placement_grows_all_conformers_and_preserves_every_old_value() {
    let result = one_hydrogen_result();
    let result_snapshot = result.clone();
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(7, vec![[1.25, -2.5]]).with_prop("source", "two")],
        conformers_3d: vec![
            Conformer3D::new(11, vec![[3.0, 4.0, 5.0]], true).with_prop("source", "three"),
            Conformer3D::new(12, vec![[-1.0, 2.0, 0.0]], false).with_prop("flat", "yes"),
        ],
        source_coordinate_dim: Some(CoordinateDimension::ThreeD),
    };
    let coordinates_snapshot = coordinates.clone();
    let output =
        add_hydrogen_coordinates(result.clone(), coordinates.clone(), false, false).unwrap();

    assert_eq!(result, result_snapshot);
    assert_eq!(coordinates, coordinates_snapshot);
    assert_eq!(output.topology, result.topology);
    assert_eq!(output.mapping, result.mapping);
    assert_eq!(output.additions, result.additions);
    assert_eq!(
        output.coordinates.source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
    assert_eq!(output.coordinates.conformers_2d[0].id(), 7);
    assert_eq!(
        output.coordinates.conformers_2d[0]
            .props()
            .get("source")
            .map(String::as_str),
        Some("two")
    );
    assert_eq!(
        output.coordinates.conformers_2d[0].coordinates(),
        &[[1.25, -2.5], [0.0, 0.0]]
    );
    assert_eq!(output.coordinates.conformers_3d[0].id(), 11);
    assert!(output.coordinates.conformers_3d[0].is_3d());
    assert_eq!(
        output.coordinates.conformers_3d[0].coordinates(),
        &[[3.0, 4.0, 5.0], [0.0, 0.0, 0.0]]
    );
    assert_eq!(output.coordinates.conformers_3d[1].id(), 12);
    assert!(!output.coordinates.conformers_3d[1].is_3d());
    assert_eq!(
        output.coordinates.conformers_3d[1]
            .props()
            .get("flat")
            .map(String::as_str),
        Some("yes")
    );
    assert_eq!(
        output.coordinates.conformers_3d[1].coordinates(),
        &[[-1.0, 2.0, 0.0], [0.0, 0.0, 0.0]]
    );
}

#[test]
fn degree_one_resets_direction_for_mixed_flat_and_true_3d_conformers() {
    let output = add_hydrogen_coordinates(
        one_hydrogen_result(),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(1, vec![[2.0, 3.0]])],
            conformers_3d: vec![
                Conformer3D::new(2, vec![[4.0, 5.0, 6.0]], true),
                Conformer3D::new(3, vec![[7.0, 8.0, 9.0]], false),
            ],
            ..Default::default()
        },
        true,
        false,
    )
    .unwrap();

    assert_eq!(
        output.coordinates.conformers_2d[0].coordinates()[1],
        [3.0, 3.0]
    );
    assert_eq!(
        output.coordinates.conformers_3d[1].coordinates()[1],
        [8.0, 8.0, 9.0]
    );
    let true_3d = output.coordinates.conformers_3d[0].coordinates()[1];
    assert_close(true_3d[0], 4.0);
    assert_close(true_3d[1], 5.0);
    assert_close(true_3d[2], 6.0 + rdkit_rb0(1) + rdkit_rb0(6));
}

#[test]
fn sequential_hydrogens_and_degree_two_hybridization_branches_match_source_geometry() {
    let sequential = add_hydrogen_coordinates(
        explicit_addition(isolated_parent(2, Hybridization::Sp3)),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(0, vec![[0.0, 0.0]])],
            ..Default::default()
        },
        true,
        false,
    )
    .unwrap();
    let rows = sequential.coordinates.conformers_2d[0].coordinates();
    assert_eq!(rows[1], [1.0, 0.0]);
    assert_close(rows[2][0].hypot(rows[2][1]), 1.0);
    assert_close(
        rows[1][0] * rows[2][0] + rows[1][1] * rows[2][1],
        109.471_f64.to_radians().cos(),
    );

    for hybridization in [Hybridization::Sp, Hybridization::Unspecified] {
        let source = star_parent(
            1,
            AtomSpec::new(Element::C).with_hybridization(hybridization),
        );
        let output = add_hydrogen_coordinates(
            explicit_addition(source),
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(0, vec![[0.0, 0.0], [1.0, 0.0]])],
                ..Default::default()
            },
            true,
            false,
        )
        .unwrap();
        assert_eq!(
            output.coordinates.conformers_2d[0].coordinates()[2],
            [-1.0, 0.0]
        );
    }

    let coincident = add_hydrogen_coordinates(
        explicit_addition(star_parent(
            1,
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
        )),
        CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(
                9,
                vec![[2.0, -1.0, 4.0], [2.0, -1.0, 4.0]],
                true,
            )],
            ..Default::default()
        },
        true,
        false,
    )
    .unwrap();
    assert_eq!(
        coincident.coordinates.conformers_3d[0].coordinates()[2],
        [2.0, -1.0, 4.0]
    );
}

#[test]
fn sp2_aromatic_double_and_conjugated_bonds_use_the_local_plane() {
    for bond_spec in [
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single).with_conjugated(true),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Aromatic).with_aromatic(true),
    ] {
        let aromatic = bond_spec.is_aromatic();
        let mut parent = AtomSpec::new(Element::C)
            .with_hybridization(Hybridization::Sp2)
            .with_explicit_hydrogens(1);
        let mut neighbor = AtomSpec::new(Element::C);
        if aromatic {
            parent = parent.with_aromatic(true);
            neighbor = neighbor.with_aromatic(true);
        }
        let source = topology(
            vec![parent, neighbor, AtomSpec::new(Element::C)],
            vec![
                bond_spec,
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
            ],
        );
        let output = add_hydrogen_coordinates(
            explicit_addition(source),
            CoordinateBlock {
                conformers_3d: vec![Conformer3D::new(
                    0,
                    vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
                    true,
                )],
                ..Default::default()
            },
            true,
            false,
        )
        .unwrap();
        let hydrogen = output.coordinates.conformers_3d[0].coordinates()[3];
        assert!(hydrogen[1] < 0.0);
        assert_close(hydrogen[2], 0.0);
        assert_close(
            distance_3d(hydrogen, [0.0, 0.0, 0.0]),
            rdkit_rb0(1) + rdkit_rb0(6),
        );
    }
}

#[test]
fn degree_three_cancellation_sp3_rotation_and_issue_678_are_exact() {
    let source = star_parent(
        2,
        AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
    );
    let cancelled = add_hydrogen_coordinates(
        explicit_addition(source.clone()),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                0,
                vec![[2.0, 3.0], [3.0, 3.0], [1.0, 3.0]],
            )],
            ..Default::default()
        },
        true,
        false,
    )
    .unwrap();
    assert_eq!(
        cancelled.coordinates.conformers_2d[0].coordinates()[3],
        [0.0, 0.0]
    );

    let coincident = add_hydrogen_coordinates(
        explicit_addition(source),
        CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(
                0,
                vec![[2.0, 3.0, 4.0], [2.0, 3.0, 4.0], [1.0, 3.0, 4.0]],
                true,
            )],
            ..Default::default()
        },
        true,
        false,
    )
    .unwrap();
    assert_eq!(
        coincident.coordinates.conformers_3d[0].coordinates()[3],
        [2.0, 3.0, 4.0]
    );

    let rotated = add_hydrogen_coordinates(
        explicit_addition(star_parent(
            2,
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
        )),
        CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(
                0,
                vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                true,
            )],
            ..Default::default()
        },
        true,
        false,
    )
    .unwrap();
    let hydrogen = rotated.coordinates.conformers_3d[0].coordinates()[3];
    assert!(hydrogen[2].abs() > 1.0e-6);
    assert_close(
        distance_3d(hydrogen, [0.0, 0.0, 0.0]),
        rdkit_rb0(1) + rdkit_rb0(6),
    );
}

fn place_degree_four(parent: AtomSpec, neighbor_rows: [[f64; 3]; 3]) -> [f64; 3] {
    let output = add_hydrogen_coordinates(
        explicit_addition(star_parent(3, parent)),
        CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(
                0,
                std::iter::once([0.0, 0.0, 0.0])
                    .chain(neighbor_rows)
                    .collect(),
                true,
            )],
            ..Default::default()
        },
        true,
        false,
    )
    .unwrap();
    output.coordinates.conformers_3d[0].coordinates()[4]
}

#[test]
fn degree_four_planar_chiral_nonplanar_bisector_default_and_degenerate_paths_are_finite() {
    let primary = place_degree_four(
        AtomSpec::new(Element::C),
        [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]],
    );
    assert!(primary[2] > 0.0);

    let fallback_rows = [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]];
    let fallback = place_degree_four(AtomSpec::new(Element::C), fallback_rows);
    assert!(fallback[2] < 0.0);
    let chiral = place_degree_four(
        AtomSpec::new(Element::C)
            .with_chiral_tag(ChiralTag::TetrahedralCw)
            .with_prop("_CIPCode", "R")
            .unwrap(),
        fallback_rows,
    );
    assert!(chiral[2] > 0.0);

    let nonplanar = place_degree_four(
        AtomSpec::new(Element::C),
        [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
    );
    assert!(
        nonplanar
            .into_iter()
            .all(|value| value > 0.0 && value.is_finite())
    );

    let bisector = add_hydrogen_coordinates(
        explicit_addition(star_parent(3, AtomSpec::new(Element::C))),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                0,
                vec![[0.0, 0.0], [0.0, -1.0], [-1.0, 0.0], [1.0, 0.0]],
            )],
            ..Default::default()
        },
        true,
        false,
    )
    .unwrap();
    let bisector_h = bisector.coordinates.conformers_2d[0].coordinates()[4];
    assert_close(bisector_h[0], 0.0);
    assert_close(bisector_h[1], 1.0);

    let colinear = place_degree_four(
        AtomSpec::new(Element::C),
        [[-1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    );
    assert_eq!(colinear, [0.0, 0.0, 0.0]);

    let default = add_hydrogen_coordinates(
        explicit_addition(star_parent(4, AtomSpec::new(Element::C))),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                0,
                vec![[5.0, 6.0], [6.0, 6.0], [5.0, 7.0], [4.0, 6.0], [5.0, 5.0]],
            )],
            ..Default::default()
        },
        true,
        false,
    )
    .unwrap();
    assert_eq!(
        default.coordinates.conformers_2d[0].coordinates()[5],
        [0.0, 0.0]
    );
}

#[test]
fn invalid_coordinate_and_mapping_inputs_return_exact_errors_without_mutating_sources() {
    let result = one_hydrogen_result();
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(3, Vec::new())],
        ..Default::default()
    };
    let result_snapshot = result.clone();
    let coordinate_snapshot = coordinates.clone();
    assert_eq!(
        add_hydrogen_coordinates(result.clone(), coordinates.clone(), false, false),
        Err(HydrogenError::InvalidCoordinates(
            CoordinateValidationError::RowCount {
                dimension: "2D",
                conformer: 3,
                rows: 0,
                atom_count: 1,
            }
        ))
    );
    assert_eq!(result, result_snapshot);
    assert_eq!(coordinates, coordinate_snapshot);

    let mut bad = one_hydrogen_result();
    bad.mapping.atoms.old_to_new.clear();
    assert_eq!(
        add_hydrogen_coordinates(bad, CoordinateBlock::default(), false, false),
        Err(HydrogenError::InvalidAdditionPlan {
            addition: None,
            reason: "new atom count does not equal old atoms plus additions",
        })
    );

    let mut bad = one_hydrogen_result();
    bad.mapping.atoms.new_to_old.pop();
    assert_eq!(
        add_hydrogen_coordinates(bad, CoordinateBlock::default(), false, false),
        Err(HydrogenError::InvalidMapping(
            MappingValidationError::Length {
                entity: "atom",
                direction: "new-to-old",
                actual: 1,
                expected: 2,
            }
        ))
    );

    let mut bad = one_hydrogen_result();
    bad.mapping.atoms.old_to_new[0] = Some(AtomId::new(2));
    assert_eq!(
        add_hydrogen_coordinates(bad, CoordinateBlock::default(), false, false),
        Err(HydrogenError::InvalidMapping(
            MappingValidationError::OutOfRange {
                entity: "atom",
                direction: "old-to-new",
                row: 0,
                mapped: 2,
                target_count: 2,
            }
        ))
    );

    let mut bad = one_hydrogen_result();
    bad.mapping.atoms.old_to_new[0] = Some(AtomId::new(1));
    assert_eq!(
        add_hydrogen_coordinates(bad, CoordinateBlock::default(), false, false),
        Err(HydrogenError::InvalidMapping(
            MappingValidationError::InverseMismatch {
                entity: "atom",
                direction: "old-to-new",
                row: 0,
                mapped: 1,
            }
        ))
    );

    let mut bad = one_hydrogen_result();
    bad.mapping.atoms.old_to_new[0] = None;
    bad.mapping.atoms.new_to_old[0] = None;
    assert_eq!(
        add_hydrogen_coordinates(bad, CoordinateBlock::default(), false, false),
        Err(HydrogenError::InvalidAdditionPlan {
            addition: None,
            reason: "mapping is not the canonical append-only mapping",
        })
    );
}

#[test]
fn every_addition_metadata_inconsistency_is_rejected_exactly() {
    let assert_plan_error =
        |result: AddHydrogensTopologyResult, addition: Option<usize>, reason: &'static str| {
            assert_eq!(
                add_hydrogen_coordinates(result, CoordinateBlock::default(), false, false),
                Err(HydrogenError::InvalidAdditionPlan { addition, reason })
            );
        };

    let mut bad = one_hydrogen_result();
    bad.topology.atoms.pop();
    bad.topology.bonds.pop();
    bad.topology.adjacency = AdjacencyList::from_topology(1, &bad.topology.bonds);
    assert_plan_error(
        bad,
        None,
        "new atom count does not equal old atoms plus additions",
    );

    let mut bad = one_hydrogen_result();
    bad.topology.bonds.pop();
    bad.topology.adjacency = AdjacencyList::from_topology(2, &bad.topology.bonds);
    assert_plan_error(
        bad,
        None,
        "new bond count does not equal old bonds plus additions",
    );

    let mut bad = one_hydrogen_result();
    bad.additions[0].atom = AtomId::new(0);
    assert_plan_error(bad, Some(0), "atom id is not the next appended row");

    let mut bad = one_hydrogen_result();
    bad.additions[0].bond = BondId::new(1);
    assert_plan_error(bad, Some(0), "bond id is not the next appended row");

    let mut bad = one_hydrogen_result();
    bad.additions[0].parent = AtomId::new(1);
    assert_plan_error(bad, Some(0), "parent is not an original atom");

    let mut bad = one_hydrogen_result();
    bad.topology.atoms[1] = Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::HE));
    assert_plan_error(bad, Some(0), "appended atom is not hydrogen");

    let mut bad = one_hydrogen_result();
    bad.topology.atoms[1] = Atom::from_spec(
        AtomId::new(1),
        AtomSpec::new(Element::H).with_implicit_hydrogen(true),
    );
    assert_plan_error(
        bad,
        Some(0),
        "hydrogen implicit marker disagrees with addition kind",
    );

    for order in [BondOrder::Double, BondOrder::Single] {
        let mut bad = one_hydrogen_result();
        let spec = if order == BondOrder::Single {
            BondSpec::new(AtomId::new(1), AtomId::new(0), order)
        } else {
            BondSpec::new(AtomId::new(0), AtomId::new(1), order)
        };
        bad.topology.bonds[0] = Bond::from_spec(BondId::new(0), spec);
        bad.topology.adjacency =
            AdjacencyList::from_topology(bad.topology.atoms.len(), &bad.topology.bonds);
        assert_plan_error(
            bad,
            Some(0),
            "appended bond is not the ordered parent-hydrogen single bond",
        );
    }
}

fn residue_plan(
    existing_hydrogens: usize,
    parent_info: AtomPdbResidueInfo,
) -> AddHydrogensTopologyResult {
    let old_atom_count = existing_hydrogens + 1;
    let old_bond_count = existing_hydrogens;
    let mut atoms = vec![Atom::from_spec(
        AtomId::new(0),
        AtomSpec::new(Element::C).with_pdb_residue_info(parent_info),
    )];
    let mut bonds = Vec::new();
    for index in 0..existing_hydrogens {
        let atom = AtomId::new(index + 1);
        atoms.push(Atom::from_spec(atom, AtomSpec::new(Element::H)));
        bonds.push(Bond::from_spec(
            BondId::new(index),
            BondSpec::new(AtomId::new(0), atom, BondOrder::Single),
        ));
    }
    let appended_atom = AtomId::new(old_atom_count);
    atoms.push(Atom::from_spec(appended_atom, AtomSpec::new(Element::H)));
    let appended_bond = BondId::new(old_bond_count);
    bonds.push(Bond::from_spec(
        appended_bond,
        BondSpec::new(AtomId::new(0), appended_atom, BondOrder::Single),
    ));
    let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap();
    AddHydrogensTopologyResult {
        topology,
        mapping: TopologyMapping::with_appended(old_atom_count, old_bond_count, 1, 1),
        additions: vec![AddedHydrogen {
            atom: appended_atom,
            bond: appended_bond,
            parent: AtomId::new(0),
            kind: AddedHydrogenKind::Explicit,
        }],
    }
}

#[test]
fn residue_assignment_preserves_existing_info_consumes_ids_wraps_names_and_uses_source_defaults() {
    let parent = AtomPdbResidueInfo::new(" C  ", -4, "LIG", 8, "Q", true)
        .with_alt_loc("A")
        .with_insertion_code("I")
        .with_occupancy(0.25)
        .with_temp_factor(12.5)
        .with_secondary_structure(3)
        .with_segment_number(4)
        .with_monomer_class("source-class");
    let existing = AtomPdbResidueInfo::new(" HX ", 90, "LIG", 8, "Q", true)
        .with_alt_loc("E")
        .with_temp_factor(7.0);
    let zero_existing = AtomPdbResidueInfo::new(" HZ ", 0, "LIG", 8, "Q", true);
    let mut plan = residue_plan(1_000, parent.clone());
    plan.topology.atoms[2].set_pdb_residue_info(Some(existing.clone()));
    plan.topology.atoms[3].set_pdb_residue_info(Some(zero_existing.clone()));
    let snapshot = plan.clone();
    let output =
        add_hydrogen_coordinates(plan.clone(), CoordinateBlock::default(), false, true).unwrap();
    assert_eq!(plan, snapshot);
    assert_eq!(output.topology.atoms[0].pdb_residue_info(), Some(&parent));
    assert_eq!(
        output.topology.atoms[1]
            .pdb_residue_info()
            .unwrap()
            .atom_name(),
        " H1 "
    );
    assert_eq!(output.topology.atoms[2].pdb_residue_info(), Some(&existing));
    assert_eq!(
        output.topology.atoms[3].pdb_residue_info(),
        Some(&zero_existing)
    );
    assert_eq!(
        output.topology.atoms[4]
            .pdb_residue_info()
            .unwrap()
            .atom_name(),
        " H4 "
    );
    assert_eq!(
        output.topology.atoms[123]
            .pdb_residue_info()
            .unwrap()
            .atom_name(),
        "3H12"
    );
    let appended = output.topology.atoms[1_001].pdb_residue_info().unwrap();
    assert_eq!(appended.atom_name(), "1H00");
    assert_eq!(appended.serial_number(), 1_088);
    assert_eq!(appended.residue_name(), "LIG");
    assert_eq!(appended.residue_number(), 8);
    assert_eq!(appended.chain_id(), "Q");
    assert!(appended.is_hetero_atom());
    assert_eq!(appended.alt_loc(), "");
    assert_eq!(appended.insertion_code(), "");
    assert_eq!(appended.occupancy(), 1.0);
    assert_eq!(appended.temp_factor(), 0.0);
    assert_eq!(appended.secondary_structure(), 0);
    assert_eq!(appended.segment_number(), 0);
    assert_eq!(appended.monomer_class(), "");
}

#[test]
fn residue_identity_resets_only_on_number_or_chain_and_disabled_flag_is_inert() {
    let mut source = topology(
        vec![
            AtomSpec::new(Element::C)
                .with_pdb_residue_info(AtomPdbResidueInfo::new(" C1 ", 1, "ALA", 1, "A", false)),
            AtomSpec::new(Element::C)
                .with_pdb_residue_info(AtomPdbResidueInfo::new(" C2 ", 2, "GLY", 1, "A", false)),
            AtomSpec::new(Element::C)
                .with_pdb_residue_info(AtomPdbResidueInfo::new(" C3 ", 3, "SER", 2, "A", false)),
            AtomSpec::new(Element::C)
                .with_pdb_residue_info(AtomPdbResidueInfo::new(" C4 ", 4, "SER", 2, "B", false)),
        ],
        Vec::new(),
    );
    for atom in &mut source.atoms {
        atom.set_explicit_hydrogens(1);
    }
    let plan = explicit_addition(source);
    let disabled =
        add_hydrogen_coordinates(plan.clone(), CoordinateBlock::default(), false, false).unwrap();
    assert!(
        disabled.topology.atoms[4..]
            .iter()
            .all(|atom| atom.pdb_residue_info().is_none())
    );
    let output = add_hydrogen_coordinates(plan, CoordinateBlock::default(), false, true).unwrap();
    let names = output.topology.atoms[4..]
        .iter()
        .map(|atom| atom.pdb_residue_info().unwrap().atom_name())
        .collect::<Vec<_>>();
    assert_eq!(names, vec![" H1 ", " H2 ", " H1 ", " H1 "]);
}

#[test]
fn composed_detached_entrypoint_grows_rows_places_coordinates_and_preserves_properties() {
    let source = isolated_parent(1, Hybridization::Sp3);
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(4, vec![[10.0, 20.0]])],
        ..Default::default()
    };
    let properties = MoleculeProperties::default();
    let source_snapshot = source.clone();
    let coordinate_snapshot = coordinates.clone();
    let property_snapshot = properties.clone();
    let output = add_hydrogens_with_params(
        source.clone(),
        coordinates.clone(),
        properties.clone(),
        &AddHsParams {
            explicit_only: true,
            add_coords: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(source, source_snapshot);
    assert_eq!(coordinates, coordinate_snapshot);
    assert_eq!(properties, property_snapshot);
    assert_eq!(output.topology.atoms.len(), 2);
    assert_eq!(
        output.coordinates.conformers_2d[0].coordinates(),
        &[[10.0, 20.0], [11.0, 20.0]]
    );
    assert_eq!(output.properties, properties);
    output.mapping.validate_for_counts(1, 2, 0, 1).unwrap();
    assert!(output.warnings.is_empty());
}
