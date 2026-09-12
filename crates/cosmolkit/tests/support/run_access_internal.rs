use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Conformer2D, CoordinateBlock, CoordinateValidationError,
    Element, MoleculeProperties, TopologyBlock, TopologyValidationError,
};

use super::*;
use crate::ops::{
    BlockAccess, CipStatePolicy, DerivedEffects, DerivedState, MappingRequirement, MoleculeOpKind,
    OperationDomain, ParityPolicy, SemanticPreconditionSet, SupportStatus, TopologyEditKind,
};

struct TestAccess;

fn spec(
    method: &'static str,
    output: MoleculeOpOutput,
    access: BlockAccess,
) -> &'static MoleculeOpSpec {
    Box::leak(Box::new(MoleculeOpSpec {
        method,
        impl_fn: "test_impl",
        output,
        result_type: "Molecule",
        domain: OperationDomain::Topology,
        kind: MoleculeOpKind::Weak,
        topology_edit: TopologyEditKind::None,
        access,
        may_mutate: access.write(),
        auto_remap: BlockSet::NONE,
        derived_effects: DerivedEffects::new(
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
        ),
        cip_state: CipStatePolicy::Preserve,
        semantic_preconditions: SemanticPreconditionSet::NONE,
        requires_mapping: MappingRequirement::None,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    }))
}

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

fn molecule() -> Molecule {
    Molecule::from_parts(
        topology(1),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(7, vec![[1.0, 2.0]])],
            ..Default::default()
        },
        MoleculeProperties::default().with_name("source"),
    )
    .expect("test molecule is valid")
}

fn denied(operation: &'static str, block: &'static str) -> OperationError {
    OperationError::AccessDenied { operation, block }
}

fn checked_out(operation: &'static str, block: &'static str) -> OperationError {
    OperationError::BlockCheckedOut { operation, block }
}

fn not_checked_out(operation: &'static str, block: &'static str) -> OperationError {
    OperationError::BlockNotCheckedOut { operation, block }
}

fn construction_error(result: Result<OpParts<'_, TestAccess>, OperationError>) -> OperationError {
    match result {
        Ok(_) => panic!("operation context construction unexpectedly succeeded"),
        Err(error) => error,
    }
}

fn operation_error<T>(result: Result<T, OperationError>) -> OperationError {
    match result {
        Ok(_) => panic!("operation unexpectedly succeeded"),
        Err(error) => error,
    }
}

#[test]
fn constructors_are_lazy_and_reject_multiple_output_before_exposure() {
    let source = molecule();
    let single = spec(
        "lazy",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
    );
    let parts = OpParts::<TestAccess>::new(&source, single).unwrap();
    assert!(std::ptr::eq(parts.source.topology(), source.topology()));
    assert!(std::ptr::eq(
        parts.source.coordinates(),
        source.coordinates()
    ));
    assert!(std::ptr::eq(parts.source.properties(), source.properties()));
    assert!(matches!(parts.topology, WorkingBlock::Shared));
    assert!(matches!(parts.coordinates, WorkingBlock::Shared));
    assert!(matches!(parts.properties, WorkingBlock::Shared));
    assert!(matches!(parts.derived_cache, WorkingBlock::Shared));

    let multiple = spec(
        "multiple",
        MoleculeOpOutput::Multiple,
        BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
    );
    assert_eq!(
        construction_error(OpParts::<TestAccess>::new(&source, multiple)),
        OperationError::OutputMismatch {
            operation: "multiple",
            expected: MoleculeOpOutput::Single,
            actual: MoleculeOpOutput::Multiple,
        }
    );

    let mut target = molecule();
    let topology_ptr = target.topology() as *const TopologyBlock;
    {
        let parts = OpParts::<TestAccess>::new_in_place(&mut target, single).unwrap();
        assert_eq!(
            parts.source.topology() as *const TopologyBlock,
            topology_ptr
        );
    }
    assert_eq!(target.properties().name(), Some("source"));
    assert_eq!(
        construction_error(OpParts::<TestAccess>::new_in_place(&mut target, multiple,)),
        OperationError::OutputMismatch {
            operation: "multiple",
            expected: MoleculeOpOutput::Single,
            actual: MoleculeOpOutput::Multiple,
        }
    );
}

#[test]
fn topology_permissions_and_lifecycle_are_exact() {
    let source = molecule();
    let none = spec(
        "topology-none",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, none).unwrap();
    assert_eq!(
        parts.read_topology_runtime().unwrap_err(),
        denied("topology-none", "topology")
    );
    assert_eq!(
        parts.checkout_topology_runtime().unwrap_err(),
        denied("topology-none", "topology")
    );
    assert_eq!(
        parts.install_topology_runtime(topology(1)).unwrap_err(),
        denied("topology-none", "topology")
    );

    let read = spec(
        "topology-read",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::TOPOLOGY, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, read).unwrap();
    assert_eq!(parts.read_topology_runtime().unwrap().atoms.len(), 1);
    assert_eq!(
        parts.checkout_topology_runtime().unwrap_err(),
        denied("topology-read", "topology")
    );

    let write = spec(
        "topology-write",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::TOPOLOGY),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, write).unwrap();
    assert_eq!(parts.read_topology_runtime().unwrap().atoms.len(), 1);
    assert_eq!(
        parts.install_topology_runtime(topology(1)).unwrap_err(),
        not_checked_out("topology-write", "topology")
    );
    let detached = parts.checkout_topology_runtime().unwrap();
    assert_eq!(detached.atoms.len(), 1);
    assert_eq!(
        parts.read_topology_runtime().unwrap_err(),
        checked_out("topology-write", "topology")
    );
    assert_eq!(
        parts.checkout_topology_runtime().unwrap_err(),
        checked_out("topology-write", "topology")
    );
    parts.install_topology_runtime(detached).unwrap();
    assert_eq!(parts.read_topology_runtime().unwrap().atoms.len(), 1);
    assert_eq!(
        parts.install_topology_runtime(topology(1)).unwrap_err(),
        not_checked_out("topology-write", "topology")
    );
}

#[test]
fn coordinate_permissions_and_lifecycle_are_exact() {
    let source = molecule();
    let none = spec(
        "coordinates-none",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, none).unwrap();
    assert_eq!(
        parts.read_coordinates_runtime().unwrap_err(),
        denied("coordinates-none", "coordinates")
    );
    assert_eq!(
        parts.checkout_coordinates_runtime().unwrap_err(),
        denied("coordinates-none", "coordinates")
    );
    assert_eq!(
        parts
            .install_coordinates_runtime(CoordinateBlock::default())
            .unwrap_err(),
        denied("coordinates-none", "coordinates")
    );

    let read = spec(
        "coordinates-read",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::COORDINATES, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, read).unwrap();
    assert_eq!(
        parts
            .read_coordinates_runtime()
            .unwrap()
            .conformers_2d
            .len(),
        1
    );
    assert_eq!(
        parts.checkout_coordinates_runtime().unwrap_err(),
        denied("coordinates-read", "coordinates")
    );

    let write = spec(
        "coordinates-write",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::COORDINATES),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, write).unwrap();
    assert_eq!(
        parts
            .read_coordinates_runtime()
            .unwrap()
            .conformers_2d
            .len(),
        1
    );
    assert_eq!(
        parts
            .install_coordinates_runtime(CoordinateBlock::default())
            .unwrap_err(),
        not_checked_out("coordinates-write", "coordinates")
    );
    let detached = parts.checkout_coordinates_runtime().unwrap();
    assert_eq!(
        parts.read_coordinates_runtime().unwrap_err(),
        checked_out("coordinates-write", "coordinates")
    );
    assert_eq!(
        parts.checkout_coordinates_runtime().unwrap_err(),
        checked_out("coordinates-write", "coordinates")
    );
    parts.install_coordinates_runtime(detached).unwrap();
    assert_eq!(
        parts
            .read_coordinates_runtime()
            .unwrap()
            .conformers_2d
            .len(),
        1
    );
    assert_eq!(
        parts
            .install_coordinates_runtime(CoordinateBlock::default())
            .unwrap_err(),
        not_checked_out("coordinates-write", "coordinates")
    );
}

#[test]
fn property_permissions_and_lifecycle_are_exact() {
    let source = molecule();
    let none = spec(
        "properties-none",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, none).unwrap();
    assert_eq!(
        parts.read_properties_runtime().unwrap_err(),
        denied("properties-none", "properties")
    );
    assert_eq!(
        parts.checkout_properties_runtime().unwrap_err(),
        denied("properties-none", "properties")
    );
    assert_eq!(
        parts
            .install_properties_runtime(MoleculeProperties::default())
            .unwrap_err(),
        denied("properties-none", "properties")
    );

    let read = spec(
        "properties-read",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::PROPERTIES, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, read).unwrap();
    assert_eq!(
        parts.read_properties_runtime().unwrap().name(),
        Some("source")
    );
    assert_eq!(
        parts.checkout_properties_runtime().unwrap_err(),
        denied("properties-read", "properties")
    );

    let write = spec(
        "properties-write",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::PROPERTIES),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, write).unwrap();
    assert_eq!(
        parts.read_properties_runtime().unwrap().name(),
        Some("source")
    );
    assert_eq!(
        parts
            .install_properties_runtime(MoleculeProperties::default())
            .unwrap_err(),
        not_checked_out("properties-write", "properties")
    );
    let detached = parts.checkout_properties_runtime().unwrap();
    assert_eq!(
        parts.read_properties_runtime().unwrap_err(),
        checked_out("properties-write", "properties")
    );
    assert_eq!(
        parts.checkout_properties_runtime().unwrap_err(),
        checked_out("properties-write", "properties")
    );
    parts.install_properties_runtime(detached).unwrap();
    assert_eq!(
        parts.read_properties_runtime().unwrap().name(),
        Some("source")
    );
    assert_eq!(
        parts
            .install_properties_runtime(MoleculeProperties::default())
            .unwrap_err(),
        not_checked_out("properties-write", "properties")
    );
}

#[test]
fn derived_cache_permissions_and_lifecycle_are_exact() {
    let source = molecule();
    let none = spec(
        "cache-none",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, none).unwrap();
    assert_eq!(
        operation_error(parts.read_derived_cache_runtime()),
        denied("cache-none", "derived_cache")
    );
    assert_eq!(
        operation_error(parts.checkout_derived_cache_runtime()),
        denied("cache-none", "derived_cache")
    );
    assert_eq!(
        parts
            .install_derived_cache_runtime(DerivedCacheBlock::default())
            .unwrap_err(),
        denied("cache-none", "derived_cache")
    );

    let read = spec(
        "cache-read",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::DERIVED_CACHE, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, read).unwrap();
    parts.read_derived_cache_runtime().unwrap();
    assert_eq!(
        operation_error(parts.checkout_derived_cache_runtime()),
        denied("cache-read", "derived_cache")
    );

    let write = spec(
        "cache-write",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::DERIVED_CACHE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, write).unwrap();
    parts.read_derived_cache_runtime().unwrap();
    assert_eq!(
        parts
            .install_derived_cache_runtime(DerivedCacheBlock::default())
            .unwrap_err(),
        not_checked_out("cache-write", "derived_cache")
    );
    let detached = parts.checkout_derived_cache_runtime().unwrap();
    assert_eq!(
        operation_error(parts.read_derived_cache_runtime()),
        checked_out("cache-write", "derived_cache")
    );
    assert_eq!(
        operation_error(parts.checkout_derived_cache_runtime()),
        checked_out("cache-write", "derived_cache")
    );
    parts.install_derived_cache_runtime(detached).unwrap();
    parts.read_derived_cache_runtime().unwrap();
    assert_eq!(
        parts
            .install_derived_cache_runtime(DerivedCacheBlock::default())
            .unwrap_err(),
        not_checked_out("cache-write", "derived_cache")
    );
}

#[test]
fn checkout_bookkeeping_is_independent_and_sources_stay_unchanged() {
    let source = molecule();
    let all = BlockSet::TOPOLOGY
        .union(BlockSet::COORDINATES)
        .union(BlockSet::PROPERTIES)
        .union(BlockSet::DERIVED_CACHE);
    let access = spec(
        "independent",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, all),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, access).unwrap();
    let topology = parts.checkout_topology_runtime().unwrap();
    let coordinates = parts.checkout_coordinates_runtime().unwrap();
    assert_eq!(
        parts.read_properties_runtime().unwrap().name(),
        Some("source")
    );
    parts.install_topology_runtime(topology).unwrap();
    assert_eq!(parts.read_topology_runtime().unwrap().atoms.len(), 1);
    assert_eq!(
        parts.read_coordinates_runtime().unwrap_err(),
        checked_out("independent", "coordinates")
    );
    parts.install_coordinates_runtime(coordinates).unwrap();
    assert_eq!(source.num_atoms(), 1);
    assert_eq!(source.properties().name(), Some("source"));

    let mut target = molecule();
    {
        let mut parts = OpParts::<TestAccess>::new_in_place(&mut target, access).unwrap();
        let properties = parts.checkout_properties_runtime().unwrap();
        parts
            .install_properties_runtime(properties.with_name("working"))
            .unwrap();
        assert_eq!(
            parts.read_properties_runtime().unwrap().name(),
            Some("working")
        );
    }
    assert_eq!(target.properties().name(), Some("source"));
}

#[test]
fn invalid_replacements_are_rejected_without_replacing_working_or_live_state() {
    let source = molecule();
    let topology_write = spec(
        "invalid-topology",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::TOPOLOGY),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, topology_write).unwrap();
    parts.checkout_topology_runtime().unwrap();
    let invalid_topology = TopologyBlock {
        atoms: vec![atom(1)],
        bonds: Vec::new(),
        adjacency: AdjacencyList::from_topology(1, &[]),
        substance_groups: Vec::new(),
        stereo_groups: Vec::new(),
    };
    assert_eq!(
        parts
            .install_topology_runtime(invalid_topology)
            .unwrap_err(),
        OperationError::InvalidTopology(TopologyValidationError::AtomIdMismatch {
            position: 0,
            id: AtomId::new(1),
        })
    );
    assert_eq!(
        parts.read_topology_runtime().unwrap_err(),
        checked_out("invalid-topology", "topology")
    );
    assert_eq!(source.num_atoms(), 1);

    let coordinate_write = spec(
        "invalid-coordinates",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::COORDINATES),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, coordinate_write).unwrap();
    parts.checkout_coordinates_runtime().unwrap();
    let invalid_coordinates = CoordinateBlock {
        conformers_2d: vec![
            Conformer2D::new(7, vec![[0.0, 0.0]]),
            Conformer2D::new(8, vec![[0.0, 0.0], [1.0, 1.0]]),
        ],
        ..Default::default()
    };
    assert_eq!(
        parts
            .install_coordinates_runtime(invalid_coordinates)
            .unwrap_err(),
        OperationError::InvalidCoordinates(CoordinateValidationError::RowCount {
            dimension: "2D",
            conformer: 8,
            rows: 2,
            atom_count: 1,
        })
    );
    assert_eq!(
        parts.read_coordinates_runtime().unwrap_err(),
        checked_out("invalid-coordinates", "coordinates")
    );
    assert_eq!(source.coordinates().conformers_2d.len(), 1);
}
