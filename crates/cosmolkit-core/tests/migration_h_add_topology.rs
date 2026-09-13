#![cfg(feature = "op-contracts-strict")]

use cosmolkit_core::{
    __migration_hydrogens::{AddedHydrogenKind, HydrogenError, add_hydrogens_topology},
    AddHsParams,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, StereoGroup, StereoGroupKind,
    SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondOrder, Element};

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

fn bond(begin: usize, end: usize) -> BondSpec {
    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single)
}

fn propane() -> TopologyBlock {
    topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        vec![bond(0, 1), bond(1, 2)],
    )
}

fn selected_parents(
    result: &cosmolkit_core::__migration_hydrogens::AddHydrogensTopologyResult,
) -> Vec<usize> {
    result
        .additions
        .iter()
        .map(|addition| addition.parent.index())
        .collect()
}

#[test]
fn defaults_match_add_hs_parameter_source() {
    let params = AddHsParams::default();
    assert!(!params.explicit_only);
    assert!(!params.add_coords);
    assert!(!params.add_residue_info);
    assert!(!params.skip_queries);
    assert_eq!(params.only_on_atoms, None);
}

#[test]
fn empty_input_returns_valid_identity_mapping_without_mutating_source() {
    let source = TopologyBlock::default();
    let snapshot = source.clone();
    let result = add_hydrogens_topology(source.clone(), &AddHsParams::default()).unwrap();
    assert_eq!(source, snapshot);
    assert_eq!(result.topology, snapshot);
    assert!(result.additions.is_empty());
    result.mapping.validate_for_counts(0, 0, 0, 0).unwrap();
}

#[test]
fn propane_adds_source_counts_in_stable_parent_atom_and_bond_order() {
    let source = propane();
    let result = add_hydrogens_topology(source.clone(), &AddHsParams::default()).unwrap();
    assert_eq!(source.atoms.len(), 3);
    assert_eq!(source.bonds.len(), 2);
    assert_eq!(result.topology.atoms.len(), 11);
    assert_eq!(result.topology.bonds.len(), 10);
    assert_eq!(selected_parents(&result), vec![0, 0, 0, 1, 1, 2, 2, 2]);
    assert!(
        result
            .additions
            .iter()
            .all(|addition| addition.kind == AddedHydrogenKind::Implicit)
    );
    for (offset, addition) in result.additions.iter().enumerate() {
        assert_eq!(addition.atom, AtomId::new(3 + offset));
        assert_eq!(addition.bond, BondId::new(2 + offset));
        let appended_bond = &result.topology.bonds[addition.bond.index()];
        assert_eq!(appended_bond.begin(), addition.parent);
        assert_eq!(appended_bond.end(), addition.atom);
        assert_eq!(appended_bond.order(), BondOrder::Single);
    }
    assert!(result.topology.validate().is_ok());
}

#[test]
fn explicit_only_converts_explicit_rows_before_a_later_full_addition() {
    let source = topology(
        vec![AtomSpec::new(Element::C).with_explicit_hydrogens(1)],
        Vec::new(),
    );
    let explicit = add_hydrogens_topology(
        source,
        &AddHsParams {
            explicit_only: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(explicit.additions.len(), 1);
    assert_eq!(explicit.additions[0].kind, AddedHydrogenKind::Explicit);
    assert_eq!(explicit.topology.atoms[0].explicit_hydrogens(), 0);
    assert!(!explicit.topology.atoms[1].implicit_hydrogen());

    let completed = add_hydrogens_topology(explicit.topology, &AddHsParams::default()).unwrap();
    assert_eq!(completed.additions.len(), 3);
    assert!(
        completed
            .additions
            .iter()
            .all(|addition| addition.kind == AddedHydrogenKind::Implicit)
    );
    assert_eq!(completed.topology.atoms.len(), 5);
}

#[test]
fn only_on_atoms_is_idempotent_ordered_and_reports_exact_bad_id() {
    let selected = add_hydrogens_topology(
        propane(),
        &AddHsParams {
            only_on_atoms: Some(vec![AtomId::new(2), AtomId::new(0), AtomId::new(2)]),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(selected_parents(&selected), vec![0, 0, 0, 2, 2, 2]);

    let empty = add_hydrogens_topology(
        propane(),
        &AddHsParams {
            only_on_atoms: Some(Vec::new()),
            ..Default::default()
        },
    )
    .unwrap();
    assert!(empty.additions.is_empty());

    assert_eq!(
        add_hydrogens_topology(
            propane(),
            &AddHsParams {
                only_on_atoms: Some(vec![AtomId::new(1), AtomId::new(7)]),
                ..Default::default()
            },
        ),
        Err(HydrogenError::OnlyOnAtomOutOfRange {
            atom: AtomId::new(7),
            atom_count: 3,
        })
    );
}

#[test]
fn skip_queries_obeys_atom_marker_for_true_and_false_modes() {
    let marked = AtomSpec::new(Element::C)
        .with_prop("_MolFileAtomQuery", "1")
        .unwrap();
    let source = topology(vec![marked], Vec::new());
    let skipped = add_hydrogens_topology(
        source.clone(),
        &AddHsParams {
            skip_queries: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(skipped.additions.is_empty());
    let included = add_hydrogens_topology(source, &AddHsParams::default()).unwrap();
    assert_eq!(included.additions.len(), 4);
}

#[test]
fn skip_queries_obeys_each_incident_bond_marker() {
    for key in ["_MolFileBondQuery", "_MolFileBondQueryComplex"] {
        let query_bond = bond(0, 1).with_prop(key, "1").unwrap();
        let source = topology(
            vec![AtomSpec::new(Element::C), AtomSpec::new(Element::C)],
            vec![query_bond],
        );
        let skipped = add_hydrogens_topology(
            source.clone(),
            &AddHsParams {
                skip_queries: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(skipped.additions.is_empty(), "marker {key}");
        let included = add_hydrogens_topology(source, &AddHsParams::default()).unwrap();
        assert_eq!(included.additions.len(), 6, "marker {key}");
    }
}

#[test]
fn explicit_and_implicit_metadata_props_and_no_implicit_are_exact() {
    let selected = AtomSpec::new(Element::C)
        .with_explicit_hydrogens(1)
        .with_no_implicit(true)
        .with_tracked_isotopic_hydrogens(vec![2, 3])
        .with_prop("ordinary", "kept")
        .unwrap()
        .with_computed_prop("computed", "drop")
        .unwrap();
    let source = topology(vec![selected], Vec::new());
    let result = add_hydrogens_topology(source, &AddHsParams::default()).unwrap();
    assert_eq!(result.additions.len(), 1);
    assert_eq!(result.additions[0].kind, AddedHydrogenKind::Explicit);
    assert!(!result.topology.atoms[1].implicit_hydrogen());
    let parent = &result.topology.atoms[0];
    assert!(parent.no_implicit());
    assert_eq!(parent.tracked_isotopic_hydrogens(), &[2, 3]);
    assert_eq!(parent.prop("ordinary"), Some("kept"));
    assert_eq!(parent.prop("computed"), None);
}

#[test]
fn aromatic_explicit_nh_preserves_aromatic_and_no_implicit_state() {
    let source = topology(
        vec![
            AtomSpec::new(Element::N)
                .with_aromatic(true)
                .with_explicit_hydrogens(1)
                .with_no_implicit(true),
        ],
        Vec::new(),
    );
    let result = add_hydrogens_topology(
        source,
        &AddHsParams {
            explicit_only: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(result.additions.len(), 1);
    assert!(result.topology.atoms[0].is_aromatic());
    assert!(result.topology.atoms[0].no_implicit());
    assert_eq!(result.topology.atoms[0].explicit_hydrogens(), 0);
}

#[test]
fn old_properties_sgroup_and_stereo_group_rows_are_preserved() {
    let atoms = vec![
        AtomSpec::new(Element::C)
            .with_explicit_hydrogens(1)
            .with_prop("atom-key", "atom-value")
            .unwrap(),
        AtomSpec::new(Element::C),
    ];
    let bonds = vec![bond(0, 1).with_prop("bond-key", "bond-value").unwrap()];
    let mut source = topology(atoms, bonds);
    source.substance_groups = vec![
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![AtomId::new(0)])
            .with_bonds(vec![BondId::new(0)]),
    ];
    source.stereo_groups = vec![StereoGroup::new(
        StereoGroupKind::Absolute,
        vec![AtomId::new(0)],
        vec![BondId::new(0)],
    )];
    source.validate().unwrap();
    let expected_sgroups = source.substance_groups.clone();
    let expected_stereo = source.stereo_groups.clone();
    let result = add_hydrogens_topology(
        source,
        &AddHsParams {
            explicit_only: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(result.topology.substance_groups, expected_sgroups);
    assert_eq!(result.topology.stereo_groups, expected_stereo);
    assert_eq!(
        result.topology.atoms[0].prop("atom-key"),
        Some("atom-value")
    );
    assert_eq!(
        result.topology.bonds[0].prop("bond-key"),
        Some("bond-value")
    );
}

#[test]
fn append_mapping_is_complete_valid_and_inverse_consistent() {
    let result = add_hydrogens_topology(propane(), &AddHsParams::default()).unwrap();
    result.mapping.validate_for_counts(3, 11, 2, 10).unwrap();
    assert_eq!(
        result.mapping.atoms.old_to_new,
        vec![
            Some(AtomId::new(0)),
            Some(AtomId::new(1)),
            Some(AtomId::new(2)),
        ]
    );
    assert_eq!(
        &result.mapping.atoms.new_to_old[..3],
        &[
            Some(AtomId::new(0)),
            Some(AtomId::new(1)),
            Some(AtomId::new(2)),
        ]
    );
    assert!(
        result.mapping.atoms.new_to_old[3..]
            .iter()
            .all(Option::is_none)
    );
    assert_eq!(
        result.mapping.bonds.old_to_new,
        vec![Some(BondId::new(0)), Some(BondId::new(1))]
    );
    assert!(
        result.mapping.bonds.new_to_old[2..]
            .iter()
            .all(Option::is_none)
    );
    assert_eq!(result.additions.len(), 8);
    assert_eq!(
        result.topology.atoms.len() - 3,
        result.topology.bonds.len() - 2
    );
}

#[test]
fn invalid_topology_returns_structured_error_without_touching_source_value() {
    let source = TopologyBlock {
        atoms: vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
        adjacency: AdjacencyList::from_topology(0, &[]),
        ..TopologyBlock::default()
    };
    let snapshot = source.clone();
    assert_eq!(
        add_hydrogens_topology(source.clone(), &AddHsParams::default()),
        Err(HydrogenError::InvalidTopology(
            TopologyValidationError::AdjacencyMismatch
        ))
    );
    assert_eq!(source, snapshot);
}
