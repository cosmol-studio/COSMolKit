#![cfg(feature = "op-contracts-strict")]

use cosmolkit_core::{
    HydrogenError, HydrogenWarning, RemoveHsParams, remove_hydrogens_impl,
    remove_hydrogens_with_params,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondDirection, BondId, BondSpec, ChiralTag, Conformer2D,
    Conformer3D, CoordinateBlock, CoordinateDimension, CoordinateValidationError,
    MoleculeProperties, SdfPropertyList, SdfPropertyListTarget, StereoGroup, StereoGroupKind,
    SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock,
};
use cosmolkit_types::{BondOrder, Element};

fn atom(index: usize) -> AtomId {
    AtomId::new(index)
}

fn bond_id(index: usize) -> BondId {
    BondId::new(index)
}

fn single(begin: usize, end: usize) -> BondSpec {
    BondSpec::new(atom(begin), atom(end), BondOrder::Single)
}

fn topology_with_state(
    atom_specs: Vec<AtomSpec>,
    bond_specs: Vec<BondSpec>,
    substance_groups: Vec<SubstanceGroup>,
    stereo_groups: Vec<StereoGroup>,
) -> TopologyBlock {
    let atoms = atom_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Atom::from_spec(atom(index), spec))
        .collect();
    let bonds = bond_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Bond::from_spec(bond_id(index), spec))
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, substance_groups, stereo_groups).unwrap()
}

fn topology(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
    topology_with_state(atom_specs, bond_specs, Vec::new(), Vec::new())
}

fn no_sanitize() -> RemoveHsParams {
    RemoveHsParams {
        sanitize: false,
        ..Default::default()
    }
}

fn carbon_hydrogen(hydrogen: AtomSpec) -> TopologyBlock {
    topology(
        vec![AtomSpec::new(Element::C), hydrogen],
        vec![single(0, 1)],
    )
}

fn overvalent_carbon(hydrogen: AtomSpec) -> TopologyBlock {
    topology(
        vec![
            AtomSpec::new(Element::C),
            hydrogen,
            AtomSpec::new(Element::F),
            AtomSpec::new(Element::F),
            AtomSpec::new(Element::F),
            AtomSpec::new(Element::F),
            AtomSpec::new(Element::F),
        ],
        vec![
            single(0, 1),
            single(0, 2),
            single(0, 3),
            single(0, 4),
            single(0, 5),
            single(0, 6),
        ],
    )
}

#[test]
fn canonical_defaults_and_remove_all_source_vector_are_complete() {
    let defaults = RemoveHsParams::default();
    assert!(!defaults.remove_degree_zero);
    assert!(!defaults.remove_higher_degrees);
    assert!(!defaults.remove_only_h_neighbors);
    assert!(!defaults.remove_isotopes);
    assert!(!defaults.remove_and_track_isotopes);
    assert!(!defaults.remove_dummy_neighbors);
    assert!(!defaults.remove_defining_bond_stereo);
    assert!(defaults.remove_with_wedged_bond);
    assert!(!defaults.remove_with_query);
    assert!(defaults.remove_mapped);
    assert!(defaults.remove_in_sgroups);
    assert!(defaults.show_warnings);
    assert!(defaults.remove_nonimplicit);
    assert!(!defaults.update_explicit_count);
    assert!(!defaults.remove_hydrides);
    assert!(!defaults.remove_nontetrahedral_neighbors);
    assert!(defaults.sanitize);

    let remove_all = RemoveHsParams {
        remove_degree_zero: true,
        remove_higher_degrees: true,
        remove_only_h_neighbors: true,
        remove_isotopes: true,
        remove_and_track_isotopes: false,
        remove_dummy_neighbors: true,
        remove_defining_bond_stereo: true,
        remove_with_wedged_bond: true,
        remove_with_query: true,
        remove_mapped: true,
        remove_in_sgroups: true,
        show_warnings: false,
        remove_nonimplicit: true,
        update_explicit_count: false,
        remove_hydrides: true,
        remove_nontetrahedral_neighbors: true,
        sanitize: false,
    };
    let protected = AtomSpec::new(Element::H)
        .with_isotope(2)
        .with_atom_map(7)
        .with_formal_charge(-1)
        .with_prop("_MolFileAtomQuery", "1")
        .unwrap();
    let output = remove_hydrogens_with_params(
        topology(vec![protected], Vec::new()),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &remove_all,
    )
    .unwrap();
    assert!(output.topology.atoms.is_empty());
    assert!(output.topology.bonds.is_empty());
    output.mapping.validate_for_counts(1, 0, 0, 0).unwrap();
    assert!(output.warnings.is_empty());
}

#[test]
fn default_and_parameterized_entrypoints_are_equivalent_and_atomic() {
    let source = carbon_hydrogen(AtomSpec::new(Element::H));
    let snapshot = source.clone();
    let short = remove_hydrogens_impl(
        source.clone(),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap();
    let explicit = remove_hydrogens_with_params(
        source.clone(),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &RemoveHsParams::default(),
    )
    .unwrap();
    assert_eq!(short, explicit);
    assert_eq!(source, snapshot);
    assert_eq!(short.topology.atoms.len(), 1);
    assert!(short.topology.bonds.is_empty());
    assert_eq!(short.mapping.atoms().old_to_new(), &[Some(atom(0)), None]);
    assert_eq!(short.mapping.atoms().new_to_old(), &[Some(atom(0))]);
    assert_eq!(short.mapping.bonds().old_to_new(), &[None]);
    assert!(short.mapping.bonds().new_to_old().is_empty());
    let final_valence = short.final_valence.as_ref().unwrap();
    assert_eq!(final_valence.explicit_valence.len(), 1);
    assert_eq!(final_valence.implicit_hydrogens.len(), 1);
}

#[test]
fn empty_candidate_uses_identity_but_clears_only_computed_properties() {
    let source = topology(
        vec![
            AtomSpec::new(Element::C)
                .with_prop("atom-user", "keep")
                .unwrap()
                .with_computed_prop("atom-cache", "drop")
                .unwrap(),
        ],
        Vec::new(),
    );
    let properties = MoleculeProperties::default()
        .with_name("named")
        .with_prop("user", "keep")
        .unwrap()
        .with_computed_prop("cache", "drop")
        .unwrap();
    let output = remove_hydrogens_with_params(
        source,
        CoordinateBlock::default(),
        properties,
        &no_sanitize(),
    )
    .unwrap();
    assert_eq!(
        output.mapping,
        cosmolkit_model::TopologyMapping::identity(1, 0)
    );
    assert_eq!(output.topology.atoms[0].prop("atom-user"), Some("keep"));
    assert_eq!(output.topology.atoms[0].prop("atom-cache"), None);
    assert_eq!(output.properties.name(), Some("named"));
    assert_eq!(output.properties.prop("user"), Some("keep"));
    assert_eq!(output.properties.prop("cache"), None);
}

#[test]
fn implicit_only_and_isotope_permissions_keep_their_source_polarity() {
    let source = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::H).with_implicit_hydrogen(true),
            AtomSpec::new(Element::H).with_isotope(2),
        ],
        vec![single(0, 1), single(0, 2), single(0, 3)],
    );
    let implicit_only = RemoveHsParams {
        remove_nonimplicit: false,
        sanitize: false,
        ..Default::default()
    };
    let output = remove_hydrogens_with_params(
        source.clone(),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &implicit_only,
    )
    .unwrap();
    assert_eq!(
        output
            .topology
            .atoms
            .iter()
            .map(|atom| (
                atom.atomic_number(),
                atom.isotope(),
                atom.implicit_hydrogen()
            ))
            .collect::<Vec<_>>(),
        vec![(6, None, false), (1, None, false), (1, Some(2), false)]
    );

    let output = remove_hydrogens_with_params(
        source,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &RemoveHsParams {
            remove_isotopes: true,
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(output.topology.atoms.len(), 1);
    assert!(output.topology.bonds.is_empty());
}

#[test]
fn two_pass_tracking_composes_original_atom_and_bond_spaces() {
    let source = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::H).with_isotope(2),
            AtomSpec::new(Element::O),
        ],
        vec![single(0, 1), single(0, 2), single(0, 3)],
    );
    let output = remove_hydrogens_with_params(
        source,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &RemoveHsParams {
            remove_and_track_isotopes: true,
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(output.topology.atoms.len(), 2);
    assert_eq!(output.topology.bonds.len(), 1);
    assert_eq!(output.topology.atoms[0].tracked_isotopic_hydrogens(), &[2]);
    assert_eq!(
        output.mapping.atoms().old_to_new(),
        &[Some(atom(0)), None, None, Some(atom(1))]
    );
    assert_eq!(
        output.mapping.atoms().new_to_old(),
        &[Some(atom(0)), Some(atom(3))]
    );
    assert_eq!(
        output.mapping.bonds().old_to_new(),
        &[None, None, Some(bond_id(0))]
    );
    assert_eq!(output.mapping.bonds().new_to_old(), &[Some(bond_id(2))]);
    output.mapping.validate_for_counts(4, 2, 3, 1).unwrap();
}

#[test]
fn coordinates_property_lists_and_row_properties_follow_one_final_mapping() {
    let source = topology(
        vec![
            AtomSpec::new(Element::C)
                .with_prop("atom-user", "c")
                .unwrap()
                .with_computed_prop("atom-cache", "drop")
                .unwrap(),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::O)
                .with_prop("atom-user", "o")
                .unwrap(),
            AtomSpec::new(Element::H),
        ],
        vec![
            single(0, 1),
            single(0, 2)
                .with_prop("bond-user", "co")
                .unwrap()
                .with_computed_prop("bond-cache", "drop")
                .unwrap(),
            single(2, 3),
        ],
    );
    let coordinates = CoordinateBlock {
        conformers_2d: vec![
            Conformer2D::new(7, vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
                .with_prop("dim", "two"),
        ],
        conformers_3d: vec![
            Conformer3D::new(
                9,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [2.0, 0.0, 0.0],
                    [3.0, 0.0, 0.0],
                ],
                true,
            )
            .with_prop("dim", "three"),
        ],
        source_coordinate_dim: Some(CoordinateDimension::ThreeD),
    };
    let properties = MoleculeProperties::default()
        .with_name("named")
        .with_sdf_data_field("raw", "preserve")
        .with_prop("user", "keep")
        .unwrap()
        .with_computed_prop("cache", "drop")
        .unwrap()
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atoms",
            vec![
                Some("c".into()),
                Some("h1".into()),
                Some("o".into()),
                Some("h2".into()),
            ],
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Bond,
            "bonds",
            vec![Some("ch".into()), Some("co".into()), Some("oh".into())],
        ));
    let output =
        remove_hydrogens_with_params(source, coordinates, properties, &no_sanitize()).unwrap();

    assert_eq!(output.coordinates.conformers_2d[0].id(), 0);
    assert_eq!(
        output.coordinates.conformers_2d[0]
            .props()
            .get("dim")
            .map(String::as_str),
        Some("two")
    );
    assert_eq!(
        output.coordinates.conformers_2d[0].coordinates(),
        &[[0.0, 0.0], [2.0, 0.0]]
    );
    assert_eq!(
        output.coordinates.conformers_3d[0].coordinates(),
        &[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
    );
    assert_eq!(
        output.coordinates.source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
    assert_eq!(output.properties.name(), Some("named"));
    assert_eq!(output.properties.prop("user"), Some("keep"));
    assert_eq!(output.properties.prop("cache"), None);
    assert_eq!(
        output.properties.sdf_data_fields(),
        &[("raw".into(), "preserve".into())]
    );
    assert_eq!(
        output.properties.sdf_property_lists()[0].values(),
        &[Some("c".into()), Some("o".into())]
    );
    assert_eq!(
        output.properties.sdf_property_lists()[1].values(),
        &[Some("co".into())]
    );
    assert_eq!(output.topology.atoms[0].prop("atom-user"), Some("c"));
    assert_eq!(output.topology.atoms[0].prop("atom-cache"), None);
    assert_eq!(output.topology.atoms[1].prop("atom-user"), Some("o"));
    assert_eq!(output.topology.bonds[0].prop("bond-user"), Some("co"));
    assert_eq!(output.topology.bonds[0].prop("bond-cache"), None);
}

#[test]
fn higher_degree_removal_remaps_sgroups_and_enhanced_stereo_without_stale_ids() {
    let groups = vec![
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![atom(0), atom(2)]),
    ];
    let stereo =
        vec![StereoGroup::new(StereoGroupKind::Or, vec![atom(0), atom(2)], Vec::new()).with_id(17)];
    let source = topology_with_state(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::O),
        ],
        vec![single(0, 1), single(1, 2)],
        groups,
        stereo,
    );
    let output = remove_hydrogens_with_params(
        source,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &RemoveHsParams {
            remove_higher_degrees: true,
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(output.topology.atoms.len(), 2);
    assert!(output.topology.bonds.is_empty());
    assert_eq!(
        output.topology.substance_groups[0].atoms(),
        &[atom(0), atom(1)]
    );
    assert_eq!(
        output.topology.stereo_groups[0].atoms(),
        &[atom(0), atom(1)]
    );
    assert_eq!(output.topology.stereo_groups[0].id(), Some(17));
    output.topology.validate().unwrap();
}

#[test]
fn sanitize_runs_only_for_outer_nonimplicit_candidate_removal() {
    let unsanitized = remove_hydrogens_with_params(
        overvalent_carbon(AtomSpec::new(Element::H)),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &no_sanitize(),
    )
    .unwrap();
    assert_eq!(unsanitized.topology.atoms.len(), 6);

    let sanitized = remove_hydrogens_with_params(
        overvalent_carbon(AtomSpec::new(Element::H)),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &RemoveHsParams::default(),
    );
    assert!(matches!(sanitized, Err(HydrogenError::Sanitize(_))));

    let implicit_only = remove_hydrogens_with_params(
        overvalent_carbon(AtomSpec::new(Element::H).with_implicit_hydrogen(true)),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &RemoveHsParams {
            remove_nonimplicit: false,
            sanitize: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(implicit_only.topology.atoms.len(), 6);

    let preliminary_only = remove_hydrogens_with_params(
        overvalent_carbon(AtomSpec::new(Element::H)),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &RemoveHsParams {
            remove_and_track_isotopes: true,
            sanitize: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(preliminary_only.topology.atoms.len(), 6);
}

#[test]
fn post_removal_chiral_explicit_h_normalization_respects_no_implicit() {
    for (no_implicit, expected) in [(false, 0), (true, 3)] {
        let center = AtomSpec::new(Element::C)
            .with_chiral_tag(ChiralTag::TetrahedralCw)
            .with_explicit_hydrogens(1)
            .with_no_implicit(no_implicit);
        let output = remove_hydrogens_with_params(
            topology(
                vec![center, AtomSpec::new(Element::H), AtomSpec::new(Element::H)],
                vec![single(0, 1), single(0, 2)],
            ),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
            &no_sanitize(),
        )
        .unwrap();
        assert_eq!(output.topology.atoms[0].explicit_hydrogens(), expected);
        // CK-VALENCE-001: intermediate chiral-H normalization remains active,
        // but sanitize=false must not manufacture a final cache assignment.
        assert!(output.final_valence.is_none());
    }
}

#[test]
fn source_warning_sites_are_typed_ordered_and_quietable() {
    let cases = [
        (
            topology(vec![AtomSpec::new(Element::H)], Vec::new()),
            RemoveHsParams {
                sanitize: false,
                ..Default::default()
            },
            HydrogenWarning::IsolatedHydrogen { hydrogen: atom(0) },
        ),
        (
            topology(
                vec![AtomSpec::new(Element::DUMMY), AtomSpec::new(Element::H)],
                vec![single(0, 1)],
            ),
            no_sanitize(),
            HydrogenWarning::DummyAtomNeighbor {
                hydrogen: atom(1),
                neighbor: atom(0),
            },
        ),
        (
            topology(
                vec![
                    AtomSpec::new(Element::PT).with_chiral_tag(ChiralTag::SquarePlanar),
                    AtomSpec::new(Element::H),
                ],
                vec![single(0, 1)],
            ),
            no_sanitize(),
            HydrogenWarning::NonTetrahedralStereoNeighbor {
                hydrogen: atom(1),
                neighbor: atom(0),
            },
        ),
        (
            topology(
                vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
                vec![single(0, 1).with_direction(BondDirection::BeginWedge)],
            ),
            RemoveHsParams {
                remove_with_wedged_bond: false,
                sanitize: false,
                ..Default::default()
            },
            HydrogenWarning::WedgedBond {
                hydrogen: atom(1),
                bond: bond_id(0),
            },
        ),
    ];
    for (source, params, expected) in cases {
        let output = remove_hydrogens_with_params(
            source.clone(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
            &params,
        )
        .unwrap();
        assert_eq!(output.topology, source);
        assert_eq!(output.warnings, vec![expected]);

        let mut quiet = params;
        quiet.show_warnings = false;
        let quiet_output = remove_hydrogens_with_params(
            source,
            CoordinateBlock::default(),
            MoleculeProperties::default(),
            &quiet,
        )
        .unwrap();
        assert!(quiet_output.warnings.is_empty());
    }
}

#[test]
fn malformed_coordinate_and_property_rows_are_exact_and_atomic() {
    let source = carbon_hydrogen(AtomSpec::new(Element::H));
    let source_snapshot = source.clone();
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(4, vec![[0.0, 0.0]])],
        ..Default::default()
    };
    let coordinates_snapshot = coordinates.clone();
    assert_eq!(
        remove_hydrogens_with_params(
            source.clone(),
            coordinates.clone(),
            MoleculeProperties::default(),
            &no_sanitize(),
        ),
        Err(HydrogenError::InvalidCoordinates(
            CoordinateValidationError::RowCount {
                dimension: "2D",
                conformer: 4,
                rows: 1,
                atom_count: 2,
            }
        ))
    );
    assert_eq!(source, source_snapshot);
    assert_eq!(coordinates, coordinates_snapshot);

    let properties = MoleculeProperties::default().with_sdf_property_list(SdfPropertyList::new(
        SdfPropertyListTarget::Bond,
        "bad",
        Vec::new(),
    ));
    let properties_snapshot = properties.clone();
    assert_eq!(
        remove_hydrogens_with_params(
            source.clone(),
            CoordinateBlock::default(),
            properties.clone(),
            &no_sanitize(),
        ),
        Err(HydrogenError::InvalidPropertyList {
            target: SdfPropertyListTarget::Bond,
            name: "bad".into(),
            expected_rows: 1,
            actual_rows: 0,
        })
    );
    assert_eq!(source, source_snapshot);
    assert_eq!(properties, properties_snapshot);
}

#[test]
fn stereo_overflow_and_sanitize_errors_propagate_without_partial_output() {
    let overflow_source = topology(
        vec![
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(u8::MAX),
            AtomSpec::new(Element::H),
        ],
        vec![single(0, 1)],
    );
    let snapshot = overflow_source.clone();
    assert_eq!(
        remove_hydrogens_with_params(
            overflow_source.clone(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
            &no_sanitize(),
        ),
        Err(HydrogenError::ExplicitHydrogenOverflow {
            atom: atom(0),
            current: u8::MAX,
        })
    );
    assert_eq!(overflow_source, snapshot);

    let sanitize_source = overvalent_carbon(AtomSpec::new(Element::H));
    let sanitize_snapshot = sanitize_source.clone();
    assert!(matches!(
        remove_hydrogens_impl(
            sanitize_source.clone(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        ),
        Err(HydrogenError::Sanitize(_))
    ));
    assert_eq!(sanitize_source, sanitize_snapshot);
}
