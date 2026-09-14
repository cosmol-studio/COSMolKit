#[path = "support/coordinate_views.rs"]
mod coordinate_views;

use std::error::Error as _;

use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, BlockSet, Bond, BondId, BondOrder, BondSpec, Conformer2D,
    Conformer3D, CoordinateBlock, Element, HydrogenError, Molecule, MoleculeOpKind,
    MoleculeOpOutput, MoleculeProperties, OperationDomain, OperationError, ParityPolicy,
    RemoveHsParams, SdfPropertyList, SdfPropertyListTarget, StateModel, StereoGroup,
    StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, SupportStatus,
    TopologyBlock, TopologyEditKind, feature_spec, operation_invariant, operation_parity,
    operation_spec, support_matrix,
};

fn atom(index: usize, spec: AtomSpec) -> Atom {
    Atom::from_spec(AtomId::new(index), spec)
}

fn bond(index: usize, begin: usize, end: usize) -> Bond {
    Bond::from_spec(
        BondId::new(index),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
    )
}

fn simple_source() -> Molecule {
    let topology = TopologyBlock::try_from_parts(
        vec![
            atom(0, AtomSpec::new(Element::C)),
            atom(1, AtomSpec::new(Element::H)),
        ],
        vec![bond(0, 0, 1)],
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    Molecule::from_parts(
        topology,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap()
}

fn no_sanitize() -> RemoveHsParams {
    RemoveHsParams {
        sanitize: false,
        ..RemoveHsParams::default()
    }
}

fn stateful_source() -> Molecule {
    let atoms = vec![
        atom(
            0,
            AtomSpec::new(Element::C)
                .with_prop("atom-note", "c")
                .unwrap()
                .with_computed_prop("_CIPCode", "R")
                .unwrap(),
        ),
        atom(1, AtomSpec::new(Element::H)),
        atom(
            2,
            AtomSpec::new(Element::O)
                .with_prop("atom-note", "o")
                .unwrap(),
        ),
        atom(3, AtomSpec::new(Element::H)),
    ];
    let bonds = vec![
        bond(0, 0, 1),
        Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single)
                .with_prop("bond-note", "co")
                .unwrap()
                .with_computed_prop("_CIPCode", "E")
                .unwrap(),
        ),
        bond(2, 2, 3),
    ];
    let substance_groups = vec![
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![AtomId::new(0), AtomId::new(2)])
            .with_bonds(vec![BondId::new(1)])
            .with_label("survivors"),
    ];
    let stereo_groups = vec![
        StereoGroup::new(
            StereoGroupKind::Or,
            vec![AtomId::new(0), AtomId::new(2)],
            vec![BondId::new(1)],
        )
        .with_id(23),
    ];
    let topology =
        TopologyBlock::try_from_parts(atoms, bonds, substance_groups, stereo_groups).unwrap();
    let coordinates = CoordinateBlock {
        conformers_2d: vec![
            Conformer2D::new(7, vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
                .with_prop("plane", "kept"),
        ],
        conformers_3d: vec![
            Conformer3D::new(
                8,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [2.0, 0.0, 0.0],
                    [3.0, 0.0, 0.0],
                ],
                true,
            )
            .with_prop("space", "kept"),
        ],
        ..CoordinateBlock::default()
    };
    let properties = MoleculeProperties::default()
        .with_name("remove-hydrogens-public")
        .with_sdf_data_field("field", "kept")
        .with_prop("ordinary", "kept")
        .unwrap()
        .with_computed_prop("_CIPComputed", "true")
        .unwrap()
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atom-values",
            vec![
                Some("c".into()),
                Some("h1".into()),
                Some("o".into()),
                Some("h2".into()),
            ],
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Bond,
            "bond-values",
            vec![Some("ch".into()), Some("co".into()), Some("oh".into())],
        ));
    Molecule::from_parts(topology, coordinates, properties).unwrap()
}

fn ring_source_with_explicit_hydrogen() -> Molecule {
    let topology = TopologyBlock::try_from_parts(
        vec![
            atom(0, AtomSpec::new(Element::C)),
            atom(1, AtomSpec::new(Element::C)),
            atom(2, AtomSpec::new(Element::C)),
            atom(3, AtomSpec::new(Element::H)),
        ],
        vec![bond(0, 0, 1), bond(1, 1, 2), bond(2, 2, 0), bond(3, 0, 3)],
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    Molecule::from_parts(
        topology,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap()
    .with_assigned_rings()
    .unwrap()
    .with_assigned_valence()
    .unwrap()
}

fn overvalent_source() -> Molecule {
    let mut atoms = vec![
        atom(0, AtomSpec::new(Element::C)),
        atom(1, AtomSpec::new(Element::H)),
    ];
    atoms.extend((2..7).map(|index| atom(index, AtomSpec::new(Element::F))));
    let bonds = (1..7)
        .enumerate()
        .map(|(index, end)| bond(index, 0, end))
        .collect();
    let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap();
    Molecule::from_parts(
        topology,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap()
}

#[test]
fn canonical_public_signatures_and_all_defaults_compile() {
    let _: fn(&Molecule) -> Result<Molecule, OperationError> = Molecule::without_hydrogens;
    let _: for<'a, 'b> fn(&'a Molecule, &'b RemoveHsParams) -> Result<Molecule, OperationError> =
        Molecule::without_hydrogens_with_params;
    let _: fn(&mut Molecule) -> Result<(), OperationError> = Molecule::remove_hydrogens_;
    let _: for<'a, 'b> fn(&'a mut Molecule, &'b RemoveHsParams) -> Result<(), OperationError> =
        Molecule::remove_hydrogens_with_params_;
    assert_eq!(
        RemoveHsParams::default(),
        RemoveHsParams {
            remove_degree_zero: false,
            remove_higher_degrees: false,
            remove_only_h_neighbors: false,
            remove_isotopes: false,
            remove_and_track_isotopes: false,
            remove_dummy_neighbors: false,
            remove_defining_bond_stereo: false,
            remove_with_wedged_bond: true,
            remove_with_query: false,
            remove_mapped: true,
            remove_in_sgroups: true,
            show_warnings: true,
            remove_nonimplicit: true,
            update_explicit_count: false,
            remove_hydrides: false,
            remove_nontetrahedral_neighbors: false,
            sanitize: true,
        }
    );
}

#[test]
fn binding_contract_has_exact_remove_hydrogen_type_and_callables() {
    let expected = [
        "types.HydrogenError",
        "types.RemoveHsParams",
        "Molecule.without_hydrogens",
        "Molecule.without_hydrogens_with_params",
        "Molecule.remove_hydrogens_",
        "Molecule.remove_hydrogens_with_params_",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| expected.contains(&row.semantic_id))
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    for row in &rows[..2] {
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::Supported);
        assert_eq!(row.parity, BindingParity::NotApplicable);
    }
    for row in &rows[2..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(row.feature, "hydrogens");
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    assert_eq!(
        rows[2].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(
        rows[3].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(rows[4].callable.unwrap().state_model, StateModel::InPlace);
    assert_eq!(rows[5].callable.unwrap().state_model, StateModel::InPlace);
}

#[test]
fn generated_registry_and_four_matrices_share_the_strong_compacting_operation() {
    let feature = feature_spec("hydrogens").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("without_hydrogens_with_params").unwrap();
    assert_eq!(spec.domain, OperationDomain::Topology);
    assert_eq!(spec.kind, MoleculeOpKind::Strong);
    assert_eq!(spec.topology_edit, TopologyEditKind::Compacting);
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.access.read(), BlockSet::NONE);
    assert_eq!(
        spec.access.write(),
        BlockSet::TOPOLOGY
            .union(BlockSet::COORDINATES)
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE)
    );
    assert_eq!(spec.may_mutate, spec.access.write());
    assert_eq!(
        spec.auto_remap,
        BlockSet::COORDINATES.union(BlockSet::PROPERTIES)
    );
    assert_eq!(format!("{:?}", spec.requires_mapping), "Required");
    assert_eq!(spec.derived_effects.recompute.bits(), 0);
    assert_eq!(spec.derived_effects.preserve.bits(), 0);
    assert_eq!(
        spec.derived_effects.invalidate.bits(),
        (1 << 0) | (1 << 1) | (1 << 3) | (1 << 4) | (1 << 6) | (1 << 7)
    );
    assert_eq!(spec.derived_effects.operation_defined.bits(), 1 << 2);
    assert_eq!(format!("{:?}", spec.cip_state), "ClearComputed");
    assert_eq!(spec.support, SupportStatus::SupportedWithRdkitParity);
    assert_eq!(spec.parity, ParityPolicy::RequiredNow);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "strong_topology_with_coordinates"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "remove_hydrogens_rdkit"
    );
    let support = support_matrix()
        .iter()
        .find(|row| {
            row.operation
                .is_some_and(|operation| std::ptr::eq(operation, spec))
        })
        .unwrap();
    assert!(std::ptr::eq(support.feature, feature));
}

#[test]
fn live_compacting_commit_remaps_rows_and_preserves_typed_references() {
    let source = stateful_source();
    let snapshot = source.clone();
    let output = source
        .without_hydrogens_with_params(&no_sanitize())
        .unwrap();

    assert_eq!(source, snapshot);
    assert_eq!(output.num_atoms(), 2);
    assert_eq!(output.num_bonds(), 1);
    assert_eq!(
        output.atoms().iter().map(Atom::id).collect::<Vec<_>>(),
        vec![AtomId::new(0), AtomId::new(1)]
    );
    assert_eq!(
        output.atoms()[0].atomic_number(),
        Element::C.atomic_number()
    );
    assert_eq!(
        output.atoms()[1].atomic_number(),
        Element::O.atomic_number()
    );
    assert_eq!(output.bonds()[0].id(), BondId::new(0));
    assert_eq!(output.bonds()[0].begin(), AtomId::new(0));
    assert_eq!(output.bonds()[0].end(), AtomId::new(1));
    assert_eq!(
        output.topology().substance_groups[0].atoms(),
        &[AtomId::new(0), AtomId::new(1)]
    );
    assert_eq!(
        output.topology().substance_groups[0].bonds(),
        &[BondId::new(0)]
    );
    assert_eq!(
        output.topology().stereo_groups[0].atoms(),
        &[AtomId::new(0), AtomId::new(1)]
    );
    assert_eq!(
        output.topology().stereo_groups[0].bonds(),
        &[BondId::new(0)]
    );
    assert_eq!(output.topology().stereo_groups[0].id(), Some(23));
    output.topology().validate().unwrap();

    assert_eq!(
        output.to_builder().coordinates().conformers_2d[0].coordinates(),
        &[[0.0, 0.0], [2.0, 0.0]]
    );
    assert_eq!(
        output.to_builder().coordinates().conformers_3d[0].coordinates(),
        &[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
    );
    assert_eq!(
        output.to_builder().coordinates().conformers_2d[0]
            .props()
            .get("plane")
            .map(String::as_str),
        Some("kept")
    );
    assert_eq!(
        output.to_builder().coordinates().conformers_3d[0]
            .props()
            .get("space")
            .map(String::as_str),
        Some("kept")
    );
    assert_eq!(output.properties().name(), Some("remove-hydrogens-public"));
    assert_eq!(output.property("ordinary"), Some("kept"));
    assert_eq!(output.property("_CIPComputed"), None);
    assert_eq!(output.atoms()[0].prop("atom-note"), Some("c"));
    assert_eq!(output.atoms()[0].prop("_CIPCode"), None);
    assert_eq!(output.atoms()[1].prop("atom-note"), Some("o"));
    assert_eq!(output.bonds()[0].prop("bond-note"), Some("co"));
    assert_eq!(output.bonds()[0].prop("_CIPCode"), None);
    assert_eq!(
        output.properties().sdf_property_lists()[0].values(),
        &[Some("c".into()), Some("o".into())]
    );
    assert_eq!(
        output.properties().sdf_property_lists()[1].values(),
        &[Some("co".into())]
    );
    assert!(!std::ptr::eq(source.topology(), output.topology()));
    coordinate_views::assert_detached_coordinates(&source, &output);
    assert!(!std::ptr::eq(source.properties(), output.properties()));
}

#[test]
fn strict_commit_installs_valence_and_invalidates_ring_state_before_recomputation() {
    let source = ring_source_with_explicit_hydrogen();
    assert!(format!("{source:?}").contains("derived_cache_is_empty: false"));
    let output = source
        .without_hydrogens_with_params(&no_sanitize())
        .unwrap();
    assert_eq!(output.num_atoms(), 3);
    assert_eq!(output.num_bonds(), 3);
    assert!(format!("{output:?}").contains("derived_cache_is_empty: false"));

    let valence_rechecked = output.with_assigned_valence().unwrap();
    assert_eq!(valence_rechecked, output);
    let rings_rechecked = output.with_assigned_rings().unwrap();
    assert_eq!(rings_rechecked.topology(), output.topology());
    assert!(std::ptr::eq(rings_rechecked.topology(), output.topology()));
    assert_eq!(
        rings_rechecked.with_assigned_rings().unwrap(),
        rings_rechecked
    );
}

#[test]
fn default_parameterized_and_inplace_forms_have_identical_single_output_order() {
    let source = simple_source();
    let short = source.without_hydrogens().unwrap();
    let explicit = source
        .without_hydrogens_with_params(&RemoveHsParams::default())
        .unwrap();
    assert_eq!(short, explicit);
    assert_eq!(short.num_atoms(), 1);
    assert_eq!(short.atoms()[0].id(), AtomId::new(0));

    let mut short_in_place = source.clone();
    short_in_place.remove_hydrogens_().unwrap();
    assert_eq!(short_in_place, short);
    let mut explicit_in_place = source.clone();
    explicit_in_place
        .remove_hydrogens_with_params_(&RemoveHsParams::default())
        .unwrap();
    assert_eq!(explicit_in_place, explicit);
    // H-remove_state is a single-output operation, so per-candidate validation is N/A;
    // the stable survivor order above is the applicable ordering contract.
}

#[test]
fn structured_sanitize_failure_is_atomic_for_value_and_inplace_wrappers() {
    let source = overvalent_source();
    let observer = source.clone();
    let error = source.without_hydrogens().unwrap_err();
    assert!(matches!(
        error,
        OperationError::Hydrogen(HydrogenError::Sanitize(_))
    ));
    assert!(error.source().is_some());
    assert_eq!(source, observer);
    assert!(std::ptr::eq(source.topology(), observer.topology()));
    coordinate_views::assert_shared_coordinates(&source, &observer);
    assert!(std::ptr::eq(source.properties(), observer.properties()));

    let mut target = source.clone();
    let target_observer = target.clone();
    let error = target.remove_hydrogens_().unwrap_err();
    assert!(matches!(
        error,
        OperationError::Hydrogen(HydrogenError::Sanitize(_))
    ));
    assert_eq!(target, source);
    assert!(std::ptr::eq(target.topology(), target_observer.topology()));
    coordinate_views::assert_shared_coordinates(&target, &target_observer);
    assert!(std::ptr::eq(
        target.properties(),
        target_observer.properties()
    ));
}
