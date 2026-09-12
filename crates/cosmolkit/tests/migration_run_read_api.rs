use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingDefault, BindingExposure, BindingItem,
    BindingKind, BindingOwner, BindingParity, BindingSupport, Bond, BondId, BondOrder, BondSpec,
    Conformer2D, Conformer3D, CoordinateBlock, CoordinateValidationError, Element, Molecule,
    MoleculeProperties, OperationError, StateModel, StereoGroup, StereoGroupKind, SubstanceGroup,
    SubstanceGroupId, SubstanceGroupKind, TopologyBlock,
};

fn atom(index: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
}

fn read_fixture() -> Molecule {
    let topology = TopologyBlock::try_from_parts(
        vec![atom(0, Element::C), atom(1, Element::O)],
        vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
        )],
        vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(1)])
                .with_bonds(vec![BondId::new(0)]),
        ],
        vec![StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(0)],
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
        .with_name("read-public")
        .with_sdf_data_field("ordered", "first")
        .with_sdf_data_field("ordered", "second")
        .with_prop("empty", "")
        .unwrap()
        .with_prop("replace", "old")
        .unwrap()
        .with_prop("replace", "new")
        .unwrap();
    Molecule::from_parts(topology, coordinates, properties).unwrap()
}

#[test]
fn all_eleven_canonical_read_signatures_compile_from_the_public_crate() {
    let _: fn(&Molecule) -> usize = Molecule::num_atoms;
    let _: fn(&Molecule) -> usize = Molecule::num_bonds;
    let _: for<'a> fn(&'a Molecule) -> &'a [Atom] = Molecule::atoms;
    let _: for<'a> fn(&'a Molecule) -> &'a [Bond] = Molecule::bonds;
    let _: for<'a> fn(&'a Molecule, AtomId) -> Option<&'a Atom> = Molecule::atom;
    let _: for<'a> fn(&'a Molecule, BondId) -> Option<&'a Bond> = Molecule::bond;
    let _: for<'a> fn(&'a Molecule) -> &'a TopologyBlock = Molecule::topology;
    let _: for<'a> fn(&'a Molecule) -> &'a CoordinateBlock = Molecule::coordinates;
    let _: for<'a> fn(&'a Molecule) -> (&'a [Conformer2D], &'a [Conformer3D]) =
        Molecule::conformers;
    let _: for<'a> fn(&'a Molecule) -> &'a MoleculeProperties = Molecule::properties;
    let _: for<'a, 'b> fn(&'a Molecule, &'b str) -> Option<&'a str> = Molecule::property;
}

#[test]
fn binding_rows_exactly_match_names_signatures_and_read_only_semantics() {
    let expected = [
        ("Molecule.num_atoms", "num_atoms", "numAtoms", "usize"),
        ("Molecule.num_bonds", "num_bonds", "numBonds", "usize"),
        ("Molecule.atoms", "atoms", "atoms", "&[crate::Atom]"),
        ("Molecule.bonds", "bonds", "bonds", "&[crate::Bond]"),
        ("Molecule.atom", "atom", "atom", "Option<&crate::Atom>"),
        ("Molecule.bond", "bond", "bond", "Option<&crate::Bond>"),
        (
            "Molecule.topology",
            "topology",
            "topology",
            "&crate::TopologyBlock",
        ),
        (
            "Molecule.coordinates",
            "coordinates",
            "coordinates",
            "&crate::CoordinateBlock",
        ),
        (
            "Molecule.conformers",
            "conformers",
            "conformers",
            "(&[crate::Conformer2D],&[crate::Conformer3D])",
        ),
        (
            "Molecule.properties",
            "properties",
            "properties",
            "&crate::MoleculeProperties",
        ),
        ("Molecule.property", "property", "property", "Option<&str>"),
    ];
    let entries = BINDING_CONTRACT
        .iter()
        .filter(|entry| expected.iter().any(|row| row.0 == entry.semantic_id))
        .collect::<Vec<_>>();
    assert_eq!(
        entries
            .iter()
            .map(|entry| entry.semantic_id)
            .collect::<Vec<_>>(),
        expected.iter().map(|row| row.0).collect::<Vec<_>>()
    );

    for (entry, (_, rust_name, javascript_name, output)) in entries.iter().zip(expected) {
        let callable = entry.callable.expect("RUN-read rows are callable");
        assert_eq!(entry.item, BindingItem::Callable);
        assert_eq!(entry.owner, BindingOwner::Molecule);
        assert_eq!(
            entry.rust_path.replace(' ', ""),
            format!("crate::Molecule::{rust_name}")
        );
        assert_eq!(entry.python_name, rust_name);
        assert_eq!(entry.javascript_name, javascript_name);
        assert_eq!(entry.feature, "runtime");
        assert_eq!(entry.exposure, BindingExposure::Public);
        assert_eq!(entry.support, BindingSupport::Supported);
        assert_eq!(entry.parity, BindingParity::NotApplicable);
        assert_eq!(callable.kind, BindingKind::Instance);
        assert_eq!(callable.output_type.replace(' ', ""), output);
        assert_eq!(callable.error_type, None);
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
    }

    for (index, name, type_name) in [
        (4, "atom_id", "crate::AtomId"),
        (5, "bond_id", "crate::BondId"),
        (10, "key", "&str"),
    ] {
        let parameter = entries[index].callable.unwrap().parameters[0];
        assert_eq!(parameter.name, name);
        assert_eq!(parameter.type_name.replace(' ', ""), type_name);
        assert_eq!(parameter.default, BindingDefault::Required);
    }
    assert!(
        entries
            .iter()
            .enumerate()
            .filter(|(index, _)| !matches!(index, 4 | 5 | 10))
            .all(|(_, entry)| entry.callable.unwrap().parameters.is_empty())
    );
}

#[test]
fn public_views_preserve_canonical_rows_typed_state_and_absence() {
    let molecule = read_fixture();
    assert_eq!(molecule.num_atoms(), 2);
    assert_eq!(molecule.num_bonds(), 1);
    assert_eq!(
        molecule.atoms().iter().map(Atom::id).collect::<Vec<_>>(),
        [AtomId::new(0), AtomId::new(1)]
    );
    assert_eq!(
        molecule.bonds().iter().map(Bond::id).collect::<Vec<_>>(),
        [BondId::new(0)]
    );
    assert!(std::ptr::eq(
        molecule.atom(AtomId::new(1)).unwrap(),
        &molecule.atoms()[1]
    ));
    assert!(std::ptr::eq(
        molecule.bond(BondId::new(0)).unwrap(),
        &molecule.bonds()[0]
    ));
    assert_eq!(molecule.atom(AtomId::new(2)), None);
    assert_eq!(molecule.bond(BondId::new(1)), None);
    assert!(std::ptr::eq(
        molecule.atoms(),
        molecule.topology().atoms.as_slice()
    ));
    assert_eq!(
        molecule.topology().substance_groups[0].atoms(),
        &[AtomId::new(1)]
    );
    assert_eq!(
        molecule.topology().stereo_groups[0].bonds(),
        &[BondId::new(0)]
    );

    let (two_d, three_d) = molecule.conformers();
    assert!(std::ptr::eq(
        two_d,
        molecule.coordinates().conformers_2d.as_slice()
    ));
    assert!(std::ptr::eq(
        three_d,
        molecule.coordinates().conformers_3d.as_slice()
    ));
    assert_eq!(two_d[0].id(), 4);
    assert_eq!(three_d[0].id(), 8);
    assert_eq!(molecule.properties().name(), Some("read-public"));
    assert_eq!(molecule.property("empty"), Some(""));
    assert_eq!(molecule.property("replace"), Some("new"));
    assert_eq!(molecule.property("missing"), None);
}

#[test]
fn invalid_detached_edit_cannot_mutate_the_live_source() {
    let source = read_fixture();
    let topology_before = source.topology().clone();
    let coordinates_before = source.coordinates().clone();
    let properties_before = source.properties().clone();
    let mut builder = source.clone().to_builder();
    builder.add_atom(AtomSpec::new(Element::N));
    let error = builder.build().unwrap_err();
    assert_eq!(
        error,
        OperationError::InvalidCoordinates(CoordinateValidationError::RowCount {
            dimension: "2D",
            conformer: 4,
            rows: 2,
            atom_count: 3,
        })
    );
    assert_eq!(source.topology(), &topology_before);
    assert_eq!(source.coordinates(), &coordinates_before);
    assert_eq!(source.properties(), &properties_before);
    assert_eq!(source.num_atoms(), 2);
}

#[test]
fn read_projection_has_no_operation_multioutput_cache_or_mutable_escape() {
    let read_methods = [
        "num_atoms",
        "num_bonds",
        "atoms",
        "bonds",
        "atom",
        "bond",
        "topology",
        "coordinates",
        "conformers",
        "properties",
        "property",
    ];
    assert!(
        cosmolkit::MOLECULE_OPS
            .iter()
            .all(|spec| !read_methods.contains(&spec.method))
    );
    assert!(
        cosmolkit::SUPPORT_MATRIX
            .iter()
            .filter_map(|row| row.operation)
            .all(|spec| !read_methods.contains(&spec.method))
    );
    assert!(
        cosmolkit::OPERATION_INVARIANT_MATRIX
            .iter()
            .all(|row| !read_methods.contains(&row.operation.method))
    );
    assert!(
        cosmolkit::PARITY_MATRIX
            .iter()
            .all(|row| !read_methods.contains(&row.operation.method))
    );

    let source = include_str!("../src/molecule.rs");
    assert_eq!(source.matches("pub struct Molecule {").count(), 1);
    for forbidden in [
        "pub fn topology_mut",
        "pub fn coordinates_mut",
        "pub fn properties_mut",
        "pub fn atom_mut",
        "pub fn bond_mut",
        "pub fn derived_cache",
        "pub fn state(",
        "pub fn read_parts",
        "pub fn install_parts",
        "pub fn emit",
    ] {
        assert!(
            !source.contains(forbidden),
            "forbidden public escape: {forbidden}"
        );
    }
}
