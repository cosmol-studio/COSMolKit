#![cfg(feature = "op-contracts-strict")]

use cosmolkit_core::{
    AddHsParams, HydrogenError, HydrogenWarning, add_hydrogens_impl, add_hydrogens_with_params,
};
use cosmolkit_model::{
    Atom, AtomId, AtomPdbResidueInfo, AtomSpec, Bond, BondId, BondOrder, BondSpec, ChiralTag,
    Conformer3D, CoordinateBlock, Element, MoleculeProperties, SdfPropertyList,
    SdfPropertyListTarget, StereoGroup, StereoGroupKind, SubstanceGroup, SubstanceGroupId,
    SubstanceGroupKind, TopologyBlock,
};

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

fn single_bond(begin: usize, end: usize) -> BondSpec {
    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single)
}

fn explicit_parent(count: u8, isotopes: Vec<u16>) -> AtomSpec {
    AtomSpec::new(Element::C)
        .with_explicit_hydrogens(count)
        .with_no_implicit(true)
        .with_tracked_isotopic_hydrogens(isotopes)
}

fn explicit_params() -> AddHsParams {
    AddHsParams {
        explicit_only: true,
        ..Default::default()
    }
}

#[test]
fn default_entrypoint_composes_topology_coordinates_properties_and_mapping() {
    let source = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::C)],
        vec![single_bond(0, 1)],
    );
    let snapshot = source.clone();
    let output = add_hydrogens_impl(
        source.clone(),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap();

    assert_eq!(source, snapshot);
    assert_eq!(output.topology.atoms.len(), 8);
    assert_eq!(output.topology.bonds.len(), 7);
    assert_eq!(output.coordinates, CoordinateBlock::default());
    assert_eq!(output.properties, MoleculeProperties::default());
    assert!(output.warnings.is_empty());
    output.mapping.validate_for_counts(2, 8, 1, 7).unwrap();
}

#[test]
fn tracked_isotopes_follow_explicit_then_implicit_addition_order_and_zero_is_absent() {
    let source = topology(
        vec![
            AtomSpec::new(Element::C)
                .with_explicit_hydrogens(1)
                .with_tracked_isotopic_hydrogens(vec![2, 0, 3, 4]),
        ],
        Vec::new(),
    );
    let output = add_hydrogens_with_params(
        source,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &AddHsParams::default(),
    )
    .unwrap();

    assert_eq!(output.topology.atoms.len(), 5);
    assert_eq!(output.topology.atoms[0].tracked_isotopic_hydrogens(), &[]);
    assert!(!output.topology.atoms[1].implicit_hydrogen());
    assert!(
        output.topology.atoms[2..]
            .iter()
            .all(Atom::implicit_hydrogen)
    );
    assert_eq!(
        output.topology.atoms[1..]
            .iter()
            .map(Atom::isotope)
            .collect::<Vec<_>>(),
        vec![Some(2), None, Some(3), Some(4)]
    );
    assert!(output.warnings.is_empty());
}

#[test]
fn short_tracked_list_labels_only_its_prefix() {
    let output = add_hydrogens_with_params(
        topology(vec![explicit_parent(2, vec![2])], Vec::new()),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &explicit_params(),
    )
    .unwrap();

    assert_eq!(output.topology.atoms[1].isotope(), Some(2));
    assert_eq!(output.topology.atoms[2].isotope(), None);
    assert!(output.warnings.is_empty());
}

#[test]
fn excess_tracked_values_are_cleared_and_reported_once_per_parent_in_order() {
    let output = add_hydrogens_with_params(
        topology(
            vec![
                explicit_parent(1, vec![2, 3, 4]),
                explicit_parent(1, vec![5, 6]),
            ],
            vec![single_bond(0, 1)],
        ),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &explicit_params(),
    )
    .unwrap();

    assert_eq!(output.topology.atoms[2].isotope(), Some(2));
    assert_eq!(output.topology.atoms[3].isotope(), Some(5));
    assert_eq!(
        output.warnings,
        vec![
            HydrogenWarning::ExtraTrackedIsotopes {
                parent: AtomId::new(0),
                count: 2,
            },
            HydrogenWarning::ExtraTrackedIsotopes {
                parent: AtomId::new(1),
                count: 1,
            },
        ]
    );
}

#[test]
fn selected_atom_additions_are_incremental_and_excluded_tracking_is_retained() {
    let source = topology(
        vec![explicit_parent(1, vec![2]), explicit_parent(1, vec![3])],
        vec![single_bond(0, 1)],
    );
    let first = add_hydrogens_with_params(
        source,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &AddHsParams {
            explicit_only: true,
            only_on_atoms: Some(vec![AtomId::new(0)]),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(first.topology.atoms[0].tracked_isotopic_hydrogens(), &[]);
    assert_eq!(first.topology.atoms[1].tracked_isotopic_hydrogens(), &[3]);
    assert_eq!(first.topology.atoms[2].isotope(), Some(2));

    let second = add_hydrogens_with_params(
        first.topology,
        first.coordinates,
        first.properties,
        &AddHsParams {
            explicit_only: true,
            only_on_atoms: Some(vec![AtomId::new(1)]),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(second.topology.atoms[2].isotope(), Some(2));
    assert_eq!(second.topology.atoms[3].isotope(), Some(3));
    assert_eq!(second.topology.atoms[1].tracked_isotopic_hydrogens(), &[]);
}

#[test]
fn skipped_query_parent_retains_explicit_count_tracking_and_computed_state() {
    let query = explicit_parent(1, vec![2])
        .with_prop("_MolFileAtomQuery", "1")
        .unwrap()
        .with_computed_prop("query-cache", "keep")
        .unwrap();
    let output = add_hydrogens_with_params(
        topology(vec![query], Vec::new()),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &AddHsParams {
            explicit_only: true,
            skip_queries: true,
            ..Default::default()
        },
    )
    .unwrap();

    assert_eq!(output.topology.atoms.len(), 1);
    assert_eq!(output.topology.atoms[0].explicit_hydrogens(), 1);
    assert_eq!(output.topology.atoms[0].tracked_isotopic_hydrogens(), &[2]);
    assert_eq!(output.topology.atoms[0].prop("query-cache"), Some("keep"));
    assert!(output.warnings.is_empty());
}

#[test]
fn molecule_properties_are_preserved_cleared_and_projected_through_append_mapping() {
    let properties = MoleculeProperties::default()
        .with_name("named")
        .with_sdf_data_field("raw", "field")
        .with_prop("ordinary", "keep")
        .unwrap()
        .with_computed_prop("computed", "drop")
        .unwrap()
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atom-values",
            vec![Some("a0".into()), None],
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Bond,
            "bond-values",
            vec![Some("b0".into())],
        ));
    let output = add_hydrogens_with_params(
        topology(
            vec![explicit_parent(1, Vec::new()), AtomSpec::new(Element::C)],
            vec![single_bond(0, 1)],
        ),
        CoordinateBlock::default(),
        properties,
        &explicit_params(),
    )
    .unwrap();

    assert_eq!(output.properties.name(), Some("named"));
    assert_eq!(
        output.properties.sdf_data_fields(),
        &[("raw".into(), "field".into())]
    );
    assert_eq!(output.properties.prop("ordinary"), Some("keep"));
    assert_eq!(output.properties.prop("computed"), None);
    assert!(output.properties.computed_prop_names().is_empty());
    let lists = output.properties.sdf_property_lists();
    assert_eq!(lists[0].values(), &[Some("a0".into()), None, None]);
    assert_eq!(lists[1].values(), &[Some("b0".into()), None]);
}

#[test]
fn every_short_or_long_property_list_is_rejected_with_exact_fields() {
    let source = topology(
        vec![explicit_parent(1, Vec::new()), AtomSpec::new(Element::C)],
        vec![single_bond(0, 1)],
    );
    for (target, rows, expected_rows) in [
        (SdfPropertyListTarget::Atom, 1, 2),
        (SdfPropertyListTarget::Atom, 3, 2),
        (SdfPropertyListTarget::Bond, 0, 1),
        (SdfPropertyListTarget::Bond, 2, 1),
    ] {
        let properties = MoleculeProperties::default()
            .with_sdf_property_list(SdfPropertyList::new(target, "bad", vec![None; rows]));
        let source_snapshot = source.clone();
        let property_snapshot = properties.clone();
        assert_eq!(
            add_hydrogens_with_params(
                source.clone(),
                CoordinateBlock::default(),
                properties.clone(),
                &explicit_params(),
            ),
            Err(HydrogenError::InvalidPropertyList {
                target,
                name: "bad".into(),
                expected_rows,
                actual_rows: rows,
            })
        );
        assert_eq!(source, source_snapshot);
        assert_eq!(properties, property_snapshot);
    }
}

#[test]
fn chiral_state_stereo_groups_sgroups_and_old_bond_properties_are_preserved() {
    let atoms = vec![
        explicit_parent(1, Vec::new())
            .with_chiral_tag(ChiralTag::TetrahedralCw)
            .with_chiral_permutation(4)
            .with_prop("atom-ordinary", "keep")
            .unwrap(),
        AtomSpec::new(Element::C),
    ];
    let bonds = vec![
        single_bond(0, 1)
            .with_prop("bond-ordinary", "keep")
            .unwrap(),
    ];
    let mut source = topology(atoms, bonds);
    source.substance_groups = vec![
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![AtomId::new(0)])
            .with_bonds(vec![BondId::new(0)]),
        SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Superatom)
            .with_atoms(vec![AtomId::new(1)])
            .with_parent(SubstanceGroupId::new(0)),
    ];
    source.stereo_groups = vec![
        StereoGroup::new(
            StereoGroupKind::Or,
            vec![AtomId::new(0)],
            vec![BondId::new(0)],
        )
        .with_id(19),
    ];
    source.validate().unwrap();
    let expected_sgroups = source.substance_groups.clone();
    let expected_stereo_groups = source.stereo_groups.clone();

    let output = add_hydrogens_with_params(
        source,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &explicit_params(),
    )
    .unwrap();
    assert_eq!(output.topology.substance_groups, expected_sgroups);
    assert_eq!(output.topology.stereo_groups, expected_stereo_groups);
    assert_eq!(
        output.topology.atoms[0].chiral_tag(),
        ChiralTag::TetrahedralCw
    );
    assert_eq!(output.topology.atoms[0].chiral_permutation(), Some(4));
    assert_eq!(output.topology.atoms[0].prop("atom-ordinary"), Some("keep"));
    assert_eq!(output.topology.bonds[0].prop("bond-ordinary"), Some("keep"));
    assert!(
        output
            .topology
            .atoms
            .iter()
            .all(|atom| !atom.no_implicit() || atom.id() == AtomId::new(0))
    );
}

#[test]
fn coordinate_and_residue_options_are_composed_without_losing_mapping() {
    let parent_info = AtomPdbResidueInfo::new(" C  ", 12, "LIG", 8, "Q", true);
    let source = topology(
        vec![explicit_parent(1, Vec::new()).with_pdb_residue_info(parent_info)],
        Vec::new(),
    );
    let coordinates = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(7, vec![[1.0, 2.0, 3.0]], true)],
        ..Default::default()
    };
    let output = add_hydrogens_with_params(
        source,
        coordinates,
        MoleculeProperties::default(),
        &AddHsParams {
            explicit_only: true,
            add_coords: true,
            add_residue_info: true,
            ..Default::default()
        },
    )
    .unwrap();

    output.mapping.validate_for_counts(1, 2, 0, 1).unwrap();
    assert_eq!(output.coordinates.conformers_3d[0].coordinates().len(), 2);
    assert_ne!(
        output.coordinates.conformers_3d[0].coordinates()[1],
        [0.0; 3]
    );
    let hydrogen_info = output.topology.atoms[1].pdb_residue_info().unwrap();
    assert_eq!(hydrogen_info.atom_name(), " H1 ");
    assert_eq!(hydrogen_info.residue_name(), "LIG");
    assert_eq!(hydrogen_info.chain_id(), "Q");
}

#[test]
fn mapping_is_identity_for_old_rows_and_none_for_each_append_in_both_directions() {
    let output = add_hydrogens_with_params(
        topology(
            vec![explicit_parent(2, Vec::new()), AtomSpec::new(Element::C)],
            vec![single_bond(0, 1)],
        ),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &explicit_params(),
    )
    .unwrap();
    output.mapping.validate_for_counts(2, 4, 1, 3).unwrap();
    assert_eq!(
        output.mapping.atoms().old_to_new(),
        &[Some(AtomId::new(0)), Some(AtomId::new(1))]
    );
    assert_eq!(
        output.mapping.atoms().new_to_old(),
        &[Some(AtomId::new(0)), Some(AtomId::new(1)), None, None]
    );
    assert_eq!(output.mapping.bonds().old_to_new(), &[Some(BondId::new(0))]);
    assert_eq!(
        output.mapping.bonds().new_to_old(),
        &[Some(BondId::new(0)), None, None]
    );
}

#[test]
fn validation_failures_return_structured_errors_without_changing_sources() {
    let source = topology(vec![explicit_parent(1, vec![2])], Vec::new());
    let source_snapshot = source.clone();
    let coordinates = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(0, vec![], true)],
        ..Default::default()
    };
    let coordinate_snapshot = coordinates.clone();
    assert!(matches!(
        add_hydrogens_with_params(
            source.clone(),
            coordinates.clone(),
            MoleculeProperties::default(),
            &explicit_params(),
        ),
        Err(HydrogenError::InvalidCoordinates(_))
    ));
    assert_eq!(source, source_snapshot);
    assert_eq!(coordinates, coordinate_snapshot);

    assert_eq!(
        add_hydrogens_with_params(
            source.clone(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
            &AddHsParams {
                only_on_atoms: Some(vec![AtomId::new(9)]),
                ..Default::default()
            },
        ),
        Err(HydrogenError::OnlyOnAtomOutOfRange {
            atom: AtomId::new(9),
            atom_count: 1,
        })
    );
    assert_eq!(source, source_snapshot);
}
