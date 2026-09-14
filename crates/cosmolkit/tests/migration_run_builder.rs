use cosmolkit::{
    BINDING_CONTRACT, BindingItem, BindingOwner, BindingParity, BindingSupport, Molecule,
    MoleculeBuilder, OperationError,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D,
    Conformer3D, CoordinateBlock, CoordinateDimension, CoordinateValidationError, Element,
    MoleculeProperties, MoleculePropertyError, SdfPropertyList, SdfPropertyListTarget, StereoGroup,
    StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock,
    TopologyEditError, TopologyValidationError,
};

fn atom(index: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
}

fn topology(elements: &[Element], bonds: &[(usize, usize, BondOrder)]) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        elements
            .iter()
            .copied()
            .enumerate()
            .map(|(index, element)| atom(index, element))
            .collect(),
        bonds
            .iter()
            .copied()
            .enumerate()
            .map(|(index, (begin, end, order))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
                )
            })
            .collect(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap()
}

#[test]
fn empty_and_complete_parts_build_through_one_checked_boundary() {
    assert_eq!(MoleculeBuilder::new().build().unwrap(), Molecule::new());

    let topology = TopologyBlock::try_from_parts(
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
    .unwrap();
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(4, vec![[0.0, 1.0], [2.0, 3.0]])],
        conformers_3d: vec![Conformer3D::new(
            8,
            vec![[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]],
            true,
        )],
        ..Default::default()
    };
    let properties = MoleculeProperties::default()
        .with_name("complete")
        .with_sdf_data_field("duplicate", "one")
        .with_sdf_data_field("duplicate", "two")
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atoms",
            vec![Some("c".into()), None],
        ));
    let molecule = MoleculeBuilder::from_parts(topology.clone(), coordinates, properties.clone())
        .build()
        .unwrap();
    assert_eq!(molecule.topology(), &topology);
    assert_eq!(molecule.properties(), &properties);
    assert_eq!(
        molecule.to_builder().coordinates().source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
}

#[test]
fn build_rejects_topology_before_coordinates_and_property_lists_exactly() {
    let invalid = TopologyBlock {
        atoms: vec![atom(1, Element::C)],
        bonds: Vec::new(),
        adjacency: AdjacencyList::from_topology(1, &[]),
        substance_groups: Vec::new(),
        stereo_groups: Vec::new(),
    };
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(3, Vec::new())],
        ..Default::default()
    };
    let properties = MoleculeProperties::default().with_sdf_property_list(SdfPropertyList::new(
        SdfPropertyListTarget::Atom,
        "bad",
        Vec::new(),
    ));
    assert_eq!(
        MoleculeBuilder::from_parts(invalid, coordinates, properties).build(),
        Err(OperationError::InvalidTopology(
            TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            }
        ))
    );

    let properties = MoleculeProperties::default().with_sdf_property_list(SdfPropertyList::new(
        SdfPropertyListTarget::Bond,
        "bond-values",
        vec![None],
    ));
    assert_eq!(
        MoleculeBuilder::from_parts(
            topology(&[Element::C], &[]),
            CoordinateBlock::default(),
            properties,
        )
        .build(),
        Err(OperationError::InvalidPropertyList {
            target: "bond",
            name: "bond-values".into(),
            values: 1,
            expected: 0,
        })
    );
}

#[test]
fn checked_topology_edits_preserve_order_and_distinguish_absence_from_bad_ids() {
    let mut builder = MoleculeBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::N));
    let a2 = builder.add_atom(AtomSpec::new(Element::O));
    let b0 = builder
        .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
        .unwrap();
    let b1 = builder
        .add_bond(BondSpec::new(a0, a2, BondOrder::Double))
        .unwrap();
    assert_eq!((b0, b1), (BondId::new(0), BondId::new(1)));
    assert_eq!(builder.degree(a0), Ok(2));
    assert_eq!(builder.neighbor_bonds(a0), Ok(vec![b0, b1]));
    assert_eq!(builder.bond_between_atoms(a1, a2), Ok(None));
    assert_eq!(builder.remove_bond_between_atoms(a1, a2), Ok(false));
    assert_eq!(
        builder.add_bond(BondSpec::new(a1, a0, BondOrder::Single)),
        Err(OperationError::InvalidTopologyEdit(
            TopologyEditError::DuplicateBond { begin: a1, end: a0 }
        ))
    );
    assert!(matches!(
        builder.add_bond(BondSpec::new(a0, AtomId::new(9), BondOrder::Single)),
        Err(OperationError::InvalidTopologyEdit(
            TopologyEditError::AtomOutOfRange { atom, atom_count: 3 }
        )) if atom == AtomId::new(9)
    ));
    assert_eq!(
        builder.add_bond(BondSpec::new(a2, a2, BondOrder::Single)),
        Err(OperationError::InvalidTopologyEdit(
            TopologyEditError::InvalidResult(TopologyValidationError::SelfLoopBond {
                bond: BondId::new(2),
                atom: a2,
            })
        ))
    );
    assert!(matches!(
        builder.degree(AtomId::new(7)),
        Err(OperationError::InvalidTopologyEdit(
            TopologyEditError::AtomOutOfRange { atom, atom_count: 3 }
        )) if atom == AtomId::new(7)
    ));
    builder.set_atom_formal_charge(a1, 1).unwrap();
    builder.set_bond_order(b1, BondOrder::Triple).unwrap();
    assert_eq!(builder.atoms()[1].formal_charge(), 1);
    assert_eq!(builder.bonds()[1].order(), BondOrder::Triple);
    assert_eq!(builder.remove_bond_between_atoms(a2, a0), Ok(true));
    assert_eq!(
        builder.bonds().iter().map(Bond::id).collect::<Vec<_>>(),
        vec![b0]
    );
}

#[test]
fn coordinate_methods_fail_closed_and_assign_dimension_local_ids() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::C));
    assert_eq!(
        builder.set_2d_coordinates(Vec::new()),
        Err(OperationError::InvalidCoordinates(
            CoordinateValidationError::RowCount {
                dimension: "2D",
                conformer: 0,
                rows: 0,
                atom_count: 1,
            }
        ))
    );
    assert!(matches!(
        builder.add_3d_conformer(vec![[0.0, f64::NAN, 0.0]]),
        Err(OperationError::InvalidCoordinates(
            CoordinateValidationError::NonFiniteCoordinate { axis: "y", .. }
        ))
    ));
    builder.set_2d_coordinates(vec![[1.0, 2.0]]).unwrap();
    assert_eq!(builder.add_2d_conformer(vec![[3.0, 4.0]]), Ok(1));
    assert_eq!(builder.add_3d_conformer(vec![[5.0, 6.0, 7.0]]), Ok(0));
    let molecule = builder.build().unwrap();
    assert_eq!(
        molecule
            .to_builder()
            .coordinates()
            .conformers_2d
            .iter()
            .map(Conformer2D::id)
            .collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert_eq!(
        molecule.to_builder().coordinates().source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
}

#[test]
fn typed_group_insertion_normalizes_sgroup_ids_and_rejects_bad_references() {
    let mut builder = MoleculeBuilder::new();
    let atom = builder.add_atom(AtomSpec::new(Element::C));
    let id = builder
        .add_substance_group(
            SubstanceGroup::new(SubstanceGroupId::new(77), SubstanceGroupKind::Data)
                .with_atoms(vec![atom]),
        )
        .unwrap();
    assert_eq!(id, SubstanceGroupId::new(0));
    assert_eq!(builder.substance_groups()[0].id(), id);
    assert_eq!(
        builder.add_stereo_group(StereoGroup::new(
            StereoGroupKind::And,
            vec![atom],
            Vec::new(),
        )),
        Ok(0)
    );
    assert!(matches!(
        builder.add_stereo_group(StereoGroup::new(
            StereoGroupKind::Or,
            vec![AtomId::new(8)],
            Vec::new(),
        )),
        Err(OperationError::InvalidTopology(
            TopologyValidationError::StereoGroupAtomOutOfRange { atom, atom_count: 1 }
        )) if atom == AtomId::new(8)
    ));
}

#[test]
fn property_configuration_preserves_order_and_reports_empty_keys() {
    assert_eq!(
        MoleculeBuilder::new().with_property(String::new(), "value".into()),
        Err(OperationError::InvalidProperty(
            MoleculePropertyError::EmptyKey
        ))
    );
    let properties = MoleculeProperties::default()
        .with_prop("ordinary", "kept")
        .unwrap()
        .with_computed_prop("computed", "cached")
        .unwrap();
    let molecule = MoleculeBuilder::new()
        .with_name("named".into())
        .with_sdf_data_field("duplicate".into(), "one".into())
        .with_sdf_data_field("duplicate".into(), "two".into())
        .with_properties(properties.clone())
        .build()
        .unwrap();
    assert_eq!(molecule.properties(), &properties);
}

#[test]
fn to_builder_edits_detached_state_without_mutating_the_source() {
    let source = MoleculeBuilder::new()
        .with_name("source".into())
        .build()
        .unwrap();
    let mut builder = source.to_builder();
    builder.add_atom(AtomSpec::new(Element::C));
    let changed = builder.with_name("changed".into()).build().unwrap();
    assert_eq!(source.num_atoms(), 0);
    assert_eq!(source.properties().name(), Some("source"));
    assert_eq!(changed.num_atoms(), 1);
    assert_eq!(changed.properties().name(), Some("changed"));
}

#[test]
fn builder_binding_and_source_guards_expose_no_bypass_or_domain_branch() {
    let builder_entries = BINDING_CONTRACT
        .iter()
        .filter(|entry| {
            entry.semantic_id == "types.MoleculeBuilder"
                || entry.semantic_id.starts_with("MoleculeBuilder.")
        })
        .collect::<Vec<_>>();
    assert_eq!(builder_entries.len(), 27);
    assert_eq!(builder_entries[0].item, BindingItem::Type);
    assert_eq!(builder_entries[0].owner, BindingOwner::Type);
    assert!(builder_entries.iter().all(|entry| {
        entry.owner == BindingOwner::Type
            && entry.support == BindingSupport::Supported
            && entry.parity == BindingParity::NotApplicable
            && entry.feature == "runtime"
    }));

    let source = include_str!("../src/molecule_builder.rs");
    assert_eq!(source.matches("pub struct MoleculeBuilder {").count(), 1);
    assert!(!source.contains("pub topology:"));
    assert!(!source.contains("pub coordinates:"));
    assert!(!source.contains("pub properties:"));
    assert!(!source.contains("topology_mut"));
    assert!(!source.contains("trust"));
    for forbidden in [
        "parse_",
        "sanitize",
        "hydrogen",
        "fingerprint",
        "descriptor",
    ] {
        assert!(
            !source.contains(forbidden),
            "forbidden domain token: {forbidden}"
        );
    }
}
