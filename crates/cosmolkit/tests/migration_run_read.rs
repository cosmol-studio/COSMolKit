use cosmolkit::{
    BINDING_CONTRACT, BindingExposure, BindingItem, BindingKind, BindingOwner, BindingParity,
    BindingSupport, Molecule, MoleculeBuilder, StateModel,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D, Conformer3D,
    CoordinateBlock, Element, MoleculeProperties, SdfPropertyList, SdfPropertyListTarget,
    StereoGroup, StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
    TopologyBlock,
};

fn atom(index: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
}

fn complete_molecule() -> Molecule {
    let topology = TopologyBlock::try_from_parts(
        vec![
            atom(0, Element::C),
            atom(1, Element::N),
            atom(2, Element::O),
        ],
        vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double),
            ),
        ],
        vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(0), AtomId::new(2)])
                .with_bonds(vec![BondId::new(1)]),
        ],
        vec![StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(1)],
            vec![BondId::new(0)],
        )],
    )
    .unwrap();
    let coordinates = CoordinateBlock {
        conformers_2d: vec![
            Conformer2D::new(4, vec![[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]]),
            Conformer2D::new(7, vec![[6.0, 7.0], [8.0, 9.0], [10.0, 11.0]]),
        ],
        conformers_3d: vec![Conformer3D::new(
            9,
            vec![[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]],
            true,
        )],
        ..Default::default()
    };
    let properties = MoleculeProperties::default()
        .with_name("read-fixture")
        .with_sdf_data_field("duplicate", "one")
        .with_sdf_data_field("duplicate", "two")
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atom-list",
            vec![Some("c".into()), None, Some("o".into())],
        ))
        .with_prop("empty", "")
        .unwrap()
        .with_prop("replace", "old")
        .unwrap()
        .with_prop("replace", "new")
        .unwrap()
        .with_computed_prop("computed", "cached")
        .unwrap();
    Molecule::from_parts(topology, coordinates, properties).unwrap()
}

#[test]
fn empty_views_are_consistent_and_all_lookups_are_absent() {
    let molecule = Molecule::new();
    assert_eq!(molecule.num_atoms(), 0);
    assert_eq!(molecule.num_bonds(), 0);
    assert!(molecule.atoms().is_empty());
    assert!(molecule.bonds().is_empty());
    assert!(molecule.atom(AtomId::new(0)).is_none());
    assert!(molecule.bond(BondId::new(0)).is_none());
    assert!(molecule.topology().substance_groups.is_empty());
    assert!(molecule.topology().stereo_groups.is_empty());
    assert_eq!(molecule.coordinates(), &CoordinateBlock::default());
    let (two_d, three_d) = molecule.conformers();
    assert!(two_d.is_empty());
    assert!(three_d.is_empty());
    assert_eq!(molecule.properties(), &MoleculeProperties::default());
    assert_eq!(molecule.property("missing"), None);
}

#[test]
fn atom_and_bond_views_preserve_rows_and_typed_lookup_bounds() {
    let molecule = complete_molecule();
    assert_eq!(molecule.num_atoms(), molecule.atoms().len());
    assert_eq!(molecule.num_bonds(), molecule.bonds().len());
    assert_eq!(
        molecule.atoms().iter().map(Atom::id).collect::<Vec<_>>(),
        vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]
    );
    assert_eq!(
        molecule.bonds().iter().map(Bond::id).collect::<Vec<_>>(),
        vec![BondId::new(0), BondId::new(1)]
    );
    for (index, row) in molecule.atoms().iter().enumerate() {
        assert!(std::ptr::eq(
            molecule.atom(AtomId::new(index)).unwrap(),
            row
        ));
    }
    for (index, row) in molecule.bonds().iter().enumerate() {
        assert!(std::ptr::eq(
            molecule.bond(BondId::new(index)).unwrap(),
            row
        ));
    }
    assert!(molecule.atom(AtomId::new(3)).is_none());
    assert!(molecule.atom(AtomId::new(usize::MAX)).is_none());
    assert!(molecule.bond(BondId::new(2)).is_none());
    assert!(molecule.bond(BondId::new(usize::MAX)).is_none());
}

#[test]
fn topology_view_is_the_same_canonical_typed_state() {
    let molecule = complete_molecule();
    assert!(std::ptr::eq(
        molecule.atoms(),
        molecule.topology().atoms.as_slice()
    ));
    assert!(std::ptr::eq(
        molecule.bonds(),
        molecule.topology().bonds.as_slice()
    ));
    assert_eq!(
        molecule.topology().substance_groups[0].atoms(),
        &[AtomId::new(0), AtomId::new(2)]
    );
    assert_eq!(
        molecule.topology().substance_groups[0].bonds(),
        &[BondId::new(1)]
    );
    assert_eq!(
        molecule.topology().stereo_groups[0].atoms(),
        &[AtomId::new(1)]
    );
    assert_eq!(
        molecule.topology().stereo_groups[0].bonds(),
        &[BondId::new(0)]
    );
}

#[test]
fn conformer_view_preserves_dimensions_order_ids_and_rows() {
    let molecule = complete_molecule();
    let (two_d, three_d) = molecule.conformers();
    assert!(std::ptr::eq(
        two_d,
        molecule.coordinates().conformers_2d.as_slice()
    ));
    assert!(std::ptr::eq(
        three_d,
        molecule.coordinates().conformers_3d.as_slice()
    ));
    assert_eq!(
        two_d.iter().map(Conformer2D::id).collect::<Vec<_>>(),
        vec![4, 7]
    );
    assert_eq!(
        three_d.iter().map(Conformer3D::id).collect::<Vec<_>>(),
        vec![9]
    );
    assert_eq!(two_d[1].coordinates()[2], [10.0, 11.0]);
    assert_eq!(three_d[0].coordinates()[2], [6.0, 7.0, 8.0]);

    let mut two_d_only = MoleculeBuilder::new();
    two_d_only.add_atom(AtomSpec::new(Element::C));
    two_d_only.add_2d_conformer(vec![[1.0, 2.0]]).unwrap();
    let two_d_only = two_d_only.build().unwrap();
    assert_eq!(two_d_only.conformers().0.len(), 1);
    assert!(two_d_only.conformers().1.is_empty());

    let mut three_d_only = MoleculeBuilder::new();
    three_d_only.add_atom(AtomSpec::new(Element::C));
    three_d_only
        .add_3d_conformer(vec![[1.0, 2.0, 3.0]])
        .unwrap();
    let three_d_only = three_d_only.build().unwrap();
    assert!(three_d_only.conformers().0.is_empty());
    assert_eq!(three_d_only.conformers().1.len(), 1);
}

#[test]
fn aggregate_and_keyed_property_views_preserve_absence_and_metadata() {
    let molecule = complete_molecule();
    assert_eq!(molecule.properties().name(), Some("read-fixture"));
    assert_eq!(
        molecule.properties().sdf_data_fields(),
        &[
            ("duplicate".into(), "one".into()),
            ("duplicate".into(), "two".into())
        ]
    );
    assert_eq!(molecule.properties().sdf_property_lists().len(), 1);
    assert_eq!(molecule.property("empty"), Some(""));
    assert_eq!(molecule.property("replace"), Some("new"));
    assert_eq!(molecule.property("computed"), Some("cached"));
    assert!(molecule.properties().is_prop_computed("computed"));
    assert_eq!(molecule.property("missing"), None);
}

#[test]
fn detached_builder_changes_cannot_mutate_read_views_of_the_source() {
    let source = complete_molecule();
    let original_atoms = source.atoms().to_vec();
    let original_bonds = source.bonds().to_vec();
    let original_coordinates = source.coordinates().clone();
    let original_properties = source.properties().clone();
    let builder = source.clone().to_builder();
    let builder = builder
        .with_property("detached".into(), "builder-only".into())
        .unwrap();
    let changed = builder.build().unwrap();
    assert_eq!(source.atoms(), original_atoms);
    assert_eq!(source.bonds(), original_bonds);
    assert_eq!(source.coordinates(), &original_coordinates);
    assert_eq!(source.properties(), &original_properties);
    assert_eq!(source.property("detached"), None);
    assert_eq!(source.num_atoms(), 3);
    assert_eq!(changed.num_atoms(), 3);
    assert_eq!(changed.property("detached"), Some("builder-only"));
}

#[test]
fn read_binding_rows_and_source_guard_expose_no_mutation_or_operation() {
    let expected_ids = [
        "Molecule.num_atoms",
        "Molecule.num_bonds",
        "Molecule.atoms",
        "Molecule.bonds",
        "Molecule.atom",
        "Molecule.bond",
        "Molecule.topology",
        "Molecule.coordinates",
        "Molecule.conformers",
        "Molecule.properties",
        "Molecule.property",
    ];
    let entries = BINDING_CONTRACT
        .iter()
        .filter(|entry| expected_ids.contains(&entry.semantic_id))
        .collect::<Vec<_>>();
    assert_eq!(
        entries
            .iter()
            .map(|entry| entry.semantic_id)
            .collect::<Vec<_>>(),
        expected_ids
    );
    assert!(entries.iter().all(|entry| {
        let callable = entry.callable.unwrap();
        entry.item == BindingItem::Callable
            && entry.owner == BindingOwner::Molecule
            && entry.exposure == BindingExposure::Public
            && entry.feature == "runtime"
            && entry.support == BindingSupport::Supported
            && entry.parity == BindingParity::NotApplicable
            && callable.kind == BindingKind::Instance
            && callable.state_model == StateModel::ReadOnly
            && callable.error_type.is_none()
            && callable.operation_semantic_id.is_none()
    }));
    assert_eq!(entries[4].callable.unwrap().parameters[0].name, "atom_id");
    assert_eq!(entries[5].callable.unwrap().parameters[0].name, "bond_id");
    assert_eq!(entries[10].callable.unwrap().parameters[0].name, "key");

    let source = include_str!("../src/molecule.rs");
    assert_eq!(source.matches("pub struct Molecule {").count(), 1);
    for forbidden in [
        "pub fn topology_mut",
        "pub fn coordinates_mut",
        "pub fn properties_mut",
        "pub fn atom_mut",
        "pub fn bond_mut",
        "pub fn state",
        "pub fn derived_cache",
        "pub fn install_parts",
        "pub fn read_parts",
    ] {
        assert!(
            !source.contains(forbidden),
            "forbidden public escape: {forbidden}"
        );
    }
}
