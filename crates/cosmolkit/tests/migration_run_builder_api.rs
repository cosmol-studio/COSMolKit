use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingDefault, BindingExposure, BindingItem,
    BindingKind, BindingOwner, BindingParity, BindingSupport, BindingTypeRole, Bond, BondId,
    BondOrder, BondSpec, Conformer2D, CoordinateBlock, CoordinateDimension,
    CoordinateValidationError, Element, Molecule, MoleculeBuilder, MoleculeProperties,
    OperationError, SdfPropertyList, SdfPropertyListTarget, StateModel, StereoGroup,
    StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock,
    TopologyValidationError,
};

fn atom(index: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
}

fn topology_with_typed_bond_dependents() -> TopologyBlock {
    TopologyBlock::try_from_parts(
        vec![atom(0, Element::C), atom(1, Element::O)],
        vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
        )],
        vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(0)])
                .with_bonds(vec![BondId::new(0)]),
        ],
        vec![StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(1)],
            vec![BondId::new(0)],
        )],
    )
    .unwrap()
}

#[test]
fn every_frozen_builder_signature_is_a_real_public_entry() {
    let _: fn() -> Molecule = Molecule::new;
    let _: fn(
        TopologyBlock,
        CoordinateBlock,
        MoleculeProperties,
    ) -> Result<Molecule, OperationError> = Molecule::from_parts;
    let _: for<'a> fn(&'a Molecule) -> MoleculeBuilder = Molecule::to_builder;

    let _: fn() -> MoleculeBuilder = MoleculeBuilder::new;
    let _: fn(TopologyBlock, CoordinateBlock, MoleculeProperties) -> MoleculeBuilder =
        MoleculeBuilder::from_parts;
    let _: fn(MoleculeBuilder) -> Result<Molecule, OperationError> = MoleculeBuilder::build;
    let _: fn(&mut MoleculeBuilder, AtomSpec) -> AtomId = MoleculeBuilder::add_atom;
    let _: fn(&mut MoleculeBuilder, BondSpec) -> Result<BondId, OperationError> =
        MoleculeBuilder::add_bond;
    let _: fn(&mut MoleculeBuilder, AtomId, i8) -> Result<(), OperationError> =
        MoleculeBuilder::set_atom_formal_charge;
    let _: fn(&mut MoleculeBuilder, BondId, BondOrder) -> Result<(), OperationError> =
        MoleculeBuilder::set_bond_order;
    let _: fn(&mut MoleculeBuilder, AtomId, AtomId) -> Result<bool, OperationError> =
        MoleculeBuilder::remove_bond_between_atoms;
    let _: fn(&MoleculeBuilder, AtomId) -> Result<usize, OperationError> = MoleculeBuilder::degree;
    let _: fn(&MoleculeBuilder, AtomId) -> Result<Vec<BondId>, OperationError> =
        MoleculeBuilder::neighbor_bonds;
    let _: fn(&MoleculeBuilder, AtomId, AtomId) -> Result<Option<BondId>, OperationError> =
        MoleculeBuilder::bond_between_atoms;
    let _: for<'a> fn(&'a MoleculeBuilder) -> &'a [Atom] = MoleculeBuilder::atoms;
    let _: for<'a> fn(&'a MoleculeBuilder) -> &'a [Bond] = MoleculeBuilder::bonds;
    let _: for<'a> fn(&'a MoleculeBuilder) -> &'a [SubstanceGroup] =
        MoleculeBuilder::substance_groups;
    let _: for<'a> fn(&'a MoleculeBuilder) -> &'a [StereoGroup] = MoleculeBuilder::stereo_groups;
    let _: for<'a> fn(&'a MoleculeBuilder) -> &'a CoordinateBlock = MoleculeBuilder::coordinates;
    let _: for<'a> fn(&'a MoleculeBuilder) -> &'a MoleculeProperties = MoleculeBuilder::properties;
    let _: fn(&mut MoleculeBuilder, Vec<[f64; 2]>) -> Result<(), OperationError> =
        MoleculeBuilder::set_2d_coordinates;
    let _: fn(&mut MoleculeBuilder, Vec<[f64; 2]>) -> Result<usize, OperationError> =
        MoleculeBuilder::add_2d_conformer;
    let _: fn(&mut MoleculeBuilder, Vec<[f64; 3]>) -> Result<usize, OperationError> =
        MoleculeBuilder::add_3d_conformer;
    let _: fn(&mut MoleculeBuilder, SubstanceGroup) -> Result<SubstanceGroupId, OperationError> =
        MoleculeBuilder::add_substance_group;
    let _: fn(&mut MoleculeBuilder, StereoGroup) -> Result<usize, OperationError> =
        MoleculeBuilder::add_stereo_group;
    let _: fn(MoleculeBuilder, String) -> MoleculeBuilder = MoleculeBuilder::with_name;
    let _: fn(MoleculeBuilder, String, String) -> Result<MoleculeBuilder, OperationError> =
        MoleculeBuilder::with_property;
    let _: fn(MoleculeBuilder, String, String) -> MoleculeBuilder =
        MoleculeBuilder::with_sdf_data_field;
    let _: fn(MoleculeBuilder, MoleculeProperties) -> MoleculeBuilder =
        MoleculeBuilder::with_properties;
}

#[test]
fn registry_exactly_matches_the_builder_receiver_defaults_feature_and_state_contract() {
    const EXPECTED_IDS: [&str; 30] = [
        "types.MoleculeBuilder",
        "Molecule.new",
        "Molecule.from_parts",
        "Molecule.to_builder",
        "MoleculeBuilder.new",
        "MoleculeBuilder.from_parts",
        "MoleculeBuilder.build",
        "MoleculeBuilder.add_atom",
        "MoleculeBuilder.add_bond",
        "MoleculeBuilder.set_atom_formal_charge",
        "MoleculeBuilder.set_bond_order",
        "MoleculeBuilder.remove_bond_between_atoms",
        "MoleculeBuilder.degree",
        "MoleculeBuilder.neighbor_bonds",
        "MoleculeBuilder.bond_between_atoms",
        "MoleculeBuilder.atoms",
        "MoleculeBuilder.bonds",
        "MoleculeBuilder.substance_groups",
        "MoleculeBuilder.stereo_groups",
        "MoleculeBuilder.coordinates",
        "MoleculeBuilder.properties",
        "MoleculeBuilder.set_2d_coordinates",
        "MoleculeBuilder.add_2d_conformer",
        "MoleculeBuilder.add_3d_conformer",
        "MoleculeBuilder.add_substance_group",
        "MoleculeBuilder.add_stereo_group",
        "MoleculeBuilder.with_name",
        "MoleculeBuilder.with_property",
        "MoleculeBuilder.with_sdf_data_field",
        "MoleculeBuilder.with_properties",
    ];
    const READ_ONLY_IDS: [&str; 8] = [
        "Molecule.to_builder",
        "MoleculeBuilder.degree",
        "MoleculeBuilder.neighbor_bonds",
        "MoleculeBuilder.bond_between_atoms",
        "MoleculeBuilder.atoms",
        "MoleculeBuilder.bonds",
        "MoleculeBuilder.substance_groups",
        "MoleculeBuilder.stereo_groups",
    ];
    const EXTRA_READ_ONLY_IDS: [&str; 2] =
        ["MoleculeBuilder.coordinates", "MoleculeBuilder.properties"];
    const IN_PLACE_IDS: [&str; 10] = [
        "MoleculeBuilder.add_atom",
        "MoleculeBuilder.add_bond",
        "MoleculeBuilder.set_atom_formal_charge",
        "MoleculeBuilder.set_bond_order",
        "MoleculeBuilder.remove_bond_between_atoms",
        "MoleculeBuilder.set_2d_coordinates",
        "MoleculeBuilder.add_2d_conformer",
        "MoleculeBuilder.add_3d_conformer",
        "MoleculeBuilder.add_substance_group",
        "MoleculeBuilder.add_stereo_group",
    ];

    let entries = BINDING_CONTRACT
        .iter()
        .filter(|entry| EXPECTED_IDS.contains(&entry.semantic_id))
        .collect::<Vec<_>>();
    assert_eq!(
        entries
            .iter()
            .map(|entry| entry.semantic_id)
            .collect::<Vec<_>>(),
        EXPECTED_IDS
    );
    assert!(entries.iter().all(|entry| {
        entry.feature == "runtime"
            && entry.exposure == BindingExposure::Public
            && entry.support == BindingSupport::Supported
            && entry.parity == BindingParity::NotApplicable
    }));

    let type_entry = entries[0];
    assert_eq!(type_entry.item, BindingItem::Type);
    assert_eq!(type_entry.owner, BindingOwner::Type);
    assert_eq!(type_entry.type_role, Some(BindingTypeRole::Value));
    assert!(type_entry.callable.is_none());

    for entry in &entries[1..] {
        let callable = entry
            .callable
            .expect("every non-type RUN-builder row is callable");
        let is_molecule_entry = entry.semantic_id.starts_with("Molecule.");
        assert_eq!(
            entry.owner,
            if is_molecule_entry {
                BindingOwner::Molecule
            } else {
                BindingOwner::Type
            }
        );
        assert_eq!(entry.item, BindingItem::Callable);
        assert!(
            callable
                .parameters
                .iter()
                .all(|parameter| parameter.default == BindingDefault::Required)
        );
        assert_eq!(callable.operation_semantic_id, None);

        let is_static = matches!(
            entry.semantic_id,
            "Molecule.new"
                | "Molecule.from_parts"
                | "MoleculeBuilder.new"
                | "MoleculeBuilder.from_parts"
        );
        assert_eq!(
            callable.kind,
            if is_static {
                BindingKind::Static
            } else {
                BindingKind::Instance
            }
        );
        let expected_state = if READ_ONLY_IDS.contains(&entry.semantic_id)
            || EXTRA_READ_ONLY_IDS.contains(&entry.semantic_id)
        {
            StateModel::ReadOnly
        } else if IN_PLACE_IDS.contains(&entry.semantic_id) {
            StateModel::InPlace
        } else {
            StateModel::ValueReturning
        };
        assert_eq!(
            callable.state_model, expected_state,
            "{}",
            entry.semantic_id
        );

        let suffix = entry.semantic_id.rsplit_once('.').unwrap().1;
        assert_eq!(entry.python_name, suffix);
    }
}

#[test]
fn detached_bond_edit_remaps_typed_dependents_and_preserves_source_value() {
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(7, vec![[0.0, 0.0], [1.0, 0.0]])],
        ..Default::default()
    };
    let properties = MoleculeProperties::default()
        .with_name("source")
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atom-values",
            vec![Some("left".into()), Some("right".into())],
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Bond,
            "bond-values",
            vec![Some("double".into())],
        ));
    let source = Molecule::from_parts(
        topology_with_typed_bond_dependents(),
        coordinates.clone(),
        properties,
    )
    .unwrap();

    let mut builder = source.to_builder();
    assert_eq!(
        builder.remove_bond_between_atoms(AtomId::new(0), AtomId::new(1)),
        Ok(true)
    );
    let changed = builder.build().unwrap();

    assert_eq!(source.num_bonds(), 1);
    assert_eq!(source.topology().substance_groups.len(), 1);
    assert_eq!(
        source.topology().stereo_groups[0].bonds(),
        &[BondId::new(0)]
    );
    assert_eq!(
        source.properties().sdf_property_lists()[1].values().len(),
        1
    );

    assert_eq!(changed.num_bonds(), 0);
    assert!(changed.topology().substance_groups.is_empty());
    assert_eq!(changed.topology().stereo_groups.len(), 1);
    assert!(changed.topology().stereo_groups[0].bonds().is_empty());
    assert_eq!(
        source.to_builder().coordinates().source_coordinate_dim,
        Some(CoordinateDimension::TwoD)
    );
    assert_eq!(
        changed.to_builder().coordinates(),
        source.to_builder().coordinates()
    );
    assert_eq!(changed.to_builder().coordinates().conformers_2d[0].id(), 7);
    assert_eq!(
        changed.to_builder().coordinates().conformers_2d[0].coordinates(),
        &[[0.0, 0.0], [1.0, 0.0]]
    );
    assert_eq!(
        changed.properties().sdf_property_lists()[0].values(),
        &[Some("left".into()), Some("right".into())]
    );
    assert!(
        changed.properties().sdf_property_lists()[1]
            .values()
            .is_empty()
    );
}

#[test]
fn invalid_construction_fails_before_authoritative_installation_with_exact_errors() {
    let invalid_topology = TopologyBlock {
        atoms: vec![atom(1, Element::C)],
        ..Default::default()
    };
    assert_eq!(
        MoleculeBuilder::from_parts(
            invalid_topology,
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .build(),
        Err(OperationError::InvalidTopology(
            TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            }
        ))
    );

    let one_atom = TopologyBlock::try_from_parts(
        vec![atom(0, Element::C)],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let invalid_coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(0, Vec::new())],
        ..Default::default()
    };
    assert_eq!(
        Molecule::from_parts(
            one_atom.clone(),
            invalid_coordinates,
            MoleculeProperties::default(),
        ),
        Err(OperationError::InvalidCoordinates(
            CoordinateValidationError::RowCount {
                dimension: "2D",
                conformer: 0,
                rows: 0,
                atom_count: 1,
            }
        ))
    );

    let invalid_properties = MoleculeProperties::default().with_sdf_property_list(
        SdfPropertyList::new(SdfPropertyListTarget::Atom, "missing-row", Vec::new()),
    );
    assert_eq!(
        MoleculeBuilder::from_parts(one_atom, CoordinateBlock::default(), invalid_properties,)
            .build(),
        Err(OperationError::InvalidPropertyList {
            target: "atom",
            name: "missing-row".into(),
            values: 0,
            expected: 1,
        })
    );

    let source = MoleculeBuilder::new()
        .with_name("unchanged".into())
        .build()
        .unwrap();
    let mut detached = source.to_builder();
    detached.add_atom(AtomSpec::new(Element::N));
    assert!(matches!(
        detached.set_2d_coordinates(Vec::new()),
        Err(OperationError::InvalidCoordinates(
            CoordinateValidationError::RowCount {
                rows: 0,
                atom_count: 1,
                ..
            }
        ))
    ));
    assert_eq!(source.num_atoms(), 0);
    assert_eq!(source.properties().name(), Some("unchanged"));
}

#[test]
fn construction_has_no_live_operation_multiple_output_cache_or_storage_bypass() {
    let entries = BINDING_CONTRACT.iter().filter(|entry| {
        entry.semantic_id == "types.MoleculeBuilder"
            || entry.semantic_id.starts_with("MoleculeBuilder.")
            || matches!(
                entry.semantic_id,
                "Molecule.new" | "Molecule.from_parts" | "Molecule.to_builder"
            )
    });
    assert!(
        entries
            .filter_map(|entry| entry.callable)
            .all(|callable| callable.operation_semantic_id.is_none())
    );

    let builder_source = include_str!("../src/molecule_builder.rs");
    assert_eq!(
        builder_source
            .matches("pub struct MoleculeBuilder {")
            .count(),
        1
    );
    for forbidden in [
        "pub topology:",
        "pub coordinates:",
        "pub properties:",
        "DerivedCacheBlock",
        "MoleculeState",
        "OpParts",
        "MultiOutput",
        "Arc<",
        "topology_mut",
        "trust",
    ] {
        assert!(
            !builder_source.contains(forbidden),
            "forbidden token: {forbidden}"
        );
    }

    let molecule_source = include_str!("../src/molecule.rs");
    let to_builder = molecule_source
        .split("pub fn to_builder")
        .nth(1)
        .unwrap()
        .split("pub fn topology")
        .next()
        .unwrap();
    assert!(!to_builder.contains("derived_cache"));
    assert!(!to_builder.contains("state.clone()"));
}
