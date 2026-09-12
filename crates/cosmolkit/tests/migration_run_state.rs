use std::mem::size_of;
use std::sync::Arc;

use cosmolkit::{
    BINDING_CONTRACT, BindingItem, BindingOwner, BindingParity, BindingSupport, BindingTypeRole,
    Molecule, OperationError,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Conformer2D, CoordinateBlock, CoordinateValidationError,
    Element, MoleculeProperties, TopologyBlock, TopologyValidationError,
};

fn atom(index: usize) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C))
}

fn topology(atom_count: usize) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        (0..atom_count).map(atom).collect(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("test topology is valid")
}

#[test]
fn empty_state_is_valid_and_has_one_pointer_sized_owner() {
    let molecule = Molecule::new();
    assert_eq!(molecule.num_atoms(), 0);
    assert_eq!(molecule.num_bonds(), 0);
    assert_eq!(molecule.topology(), &TopologyBlock::default());
    assert_eq!(molecule.coordinates(), &CoordinateBlock::default());
    assert_eq!(molecule.properties(), &MoleculeProperties::default());
    assert_eq!(size_of::<Molecule>(), size_of::<Arc<()>>());
    assert!(format!("{molecule:?}").contains("derived_cache_is_empty: true"));
}

#[test]
fn clone_shares_state_without_exposing_aliasing() {
    let source = Molecule::from_parts(
        topology(2),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(7, vec![[0.0, 1.0], [2.0, 3.0]])],
            ..Default::default()
        },
        MoleculeProperties::default().with_name("source"),
    )
    .unwrap();
    let cloned = source.clone();

    assert_eq!(cloned, source);
    assert!(std::ptr::eq(cloned.topology(), source.topology()));
    assert!(std::ptr::eq(cloned.coordinates(), source.coordinates()));
    assert!(std::ptr::eq(cloned.properties(), source.properties()));
    drop(source);
    assert_eq!(cloned.num_atoms(), 2);
    assert_eq!(cloned.properties().name(), Some("source"));
}

#[test]
fn value_equality_does_not_depend_on_arc_identity() {
    let left = Molecule::from_parts(
        topology(1),
        CoordinateBlock::default(),
        MoleculeProperties::default().with_name("same"),
    )
    .unwrap();
    let right = Molecule::from_parts(
        topology(1),
        CoordinateBlock::default(),
        MoleculeProperties::default().with_name("same"),
    )
    .unwrap();

    assert_eq!(left, right);
    assert!(!std::ptr::eq(left.topology(), right.topology()));
    assert_eq!(format!("{left:?}"), format!("{right:?}"));
}

#[test]
fn invalid_topology_is_rejected_with_exact_fields() {
    let invalid = TopologyBlock {
        atoms: vec![atom(1)],
        bonds: Vec::new(),
        adjacency: AdjacencyList::from_topology(1, &[]),
        substance_groups: Vec::new(),
        stereo_groups: Vec::new(),
    };

    assert_eq!(
        Molecule::from_parts(
            invalid,
            CoordinateBlock::default(),
            MoleculeProperties::default()
        ),
        Err(OperationError::InvalidTopology(
            TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            }
        ))
    );
}

#[test]
fn invalid_coordinate_rows_are_rejected_with_exact_fields() {
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(9, vec![[0.0, 0.0]])],
        ..Default::default()
    };

    assert_eq!(
        Molecule::from_parts(topology(2), coordinates, MoleculeProperties::default()),
        Err(OperationError::InvalidCoordinates(
            CoordinateValidationError::RowCount {
                dimension: "2D",
                conformer: 9,
                rows: 1,
                atom_count: 2,
            }
        ))
    );
}

#[test]
fn molecule_binding_is_the_only_run_state_public_entry() {
    let entries = BINDING_CONTRACT
        .iter()
        .filter(|entry| entry.semantic_id == "types.Molecule")
        .collect::<Vec<_>>();
    assert_eq!(entries.len(), 1);
    let entry = entries[0];
    assert_eq!(entry.item, BindingItem::Type);
    assert_eq!(entry.owner, BindingOwner::Type);
    assert_eq!(entry.python_name, "Molecule");
    assert_eq!(entry.javascript_name, "Molecule");
    assert_eq!(entry.feature, "runtime");
    assert_eq!(entry.support, BindingSupport::Supported);
    assert_eq!(entry.parity, BindingParity::NotApplicable);
    assert_eq!(entry.type_role, Some(BindingTypeRole::Value));
    assert_eq!(entry.callable, None);
    assert!(
        !BINDING_CONTRACT
            .iter()
            .any(|entry| entry.semantic_id.contains("MoleculeState")
                || entry.semantic_id.contains("DerivedCache"))
    );
}

#[test]
fn source_shape_keeps_state_cache_and_cow_install_private() {
    let source = include_str!("../src/molecule.rs");
    let context = include_str!("../src/ops/context.rs");
    assert_eq!(source.matches("pub struct Molecule {").count(), 1);
    assert!(source.contains("struct MoleculeState {"));
    assert!(source.contains("state: Arc<MoleculeState>"));
    assert!(source.contains("pub(crate) struct DerivedCacheBlock {"));
    for shared_block in [
        "topology: Arc<TopologyBlock>",
        "coordinates: Arc<CoordinateBlock>",
        "properties: Arc<MoleculeProperties>",
        "derived_cache: Arc<DerivedCacheBlock>",
    ] {
        assert!(source.contains(shared_block));
    }
    assert!(context.contains("let replacement = self.finish()?;"));
    assert!(context.contains("*target = replacement;"));
    assert!(!source.contains("pub struct MoleculeState"));
    assert!(!source.contains("pub state:"));
    assert!(!source.contains("pub derived_cache:"));
    assert!(!source.contains("pub fn topology_mut"));
    assert!(!source.contains("pub fn coordinates_mut"));
    assert!(!source.contains("pub fn properties_mut"));
}
